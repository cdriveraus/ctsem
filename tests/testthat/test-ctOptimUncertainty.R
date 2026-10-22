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

test_that("the normal check's bars are the coverage band's", {
  bars <- ctsem:::.ctNormalResidualBars()
  # A true sd 1.19 times the normal's leaves its 95% interval covering 90%;
  # 0.76 times, 99%; and a skew of 0.45 moves one tail to 5%.
  expect_equal(2 * pnorm(-qnorm(0.975) / bars$narrow), 0.10)
  expect_equal(2 * pnorm(-qnorm(0.975) / bars$wide), 0.01)
  drop <- qnorm(0.975)^2 / 2 - bars$skew * qnorm(0.975)^3 / 6
  expect_equal(pnorm(-sqrt(2 * drop)), 0.05)
})

# Log densities with gradients, and their exact moments on a grid.
.checkedTarget <- function(lp, grad, mode, grid) {
  lpg <- function(x) {
    out <- lp(matrix(x, 1))
    attr(out, 'gradient') <- grad(x)
    out
  }
  v <- lp(grid)
  w <- exp(v - max(v))
  w <- w / sum(w)
  mean <- colSums(grid * w)
  list(lpg = lpg, mode = mode, mean = unname(mean),
    sd = unname(sqrt(colSums(sweep(grid, 2, mean)^2 * w))))
}
.checkedRun <- function(target, uncertainty, hessian = NULL) {
  out <- suppressMessages(ctsem:::ctOptimComputeUncertainty(target$mode,
    list(nsubjects = 1), sm = NULL, lpgFunc = target$lpg,
    uncertainty = uncertainty, finishsamples = 1000, verbose = 0,
    hessian = hessian))
  drawn <- suppressWarnings(suppressMessages(ctsem:::.ctOptimDrawSamples(out,
    ctsem:::.ctResolveDraws(uncertainty, 'auto'), list(), target$mode, 1000,
    lpg = ctsem:::.ctImisPointGradbatch(target$lpg))))
  c(drawn, list(check = out$details$normalcheck))
}

test_that("checked leaves the normal of a Gaussian posterior alone", {
  set.seed(3)
  S <- crossprod(matrix(rnorm(36), 6)) / 6 + diag(6) / 10
  P <- solve(S)
  lpg <- function(x) {
    out <- -0.5 * sum(x * (P %*% x))
    attr(out, 'gradient') <- -as.numeric(P %*% x)
    out
  }
  out <- ctsem:::ctOptimComputeUncertainty(rep(0, 6), list(nsubjects = 1),
    sm = NULL, lpgFunc = lpg, uncertainty = 'checked', finishsamples = 20,
    verbose = 0)
  hess <- ctsem:::ctOptimComputeUncertainty(rep(0, 6), list(nsubjects = 1),
    sm = NULL, lpgFunc = lpg, uncertainty = 'hessian', finishsamples = 20,
    verbose = 0)
  expect_equal(ncol(out$details$normalcheck$directions), 0L)
  # Two batches of 8 that name nothing: the check stops there.
  expect_equal(out$details$normalcheck$draws, 16L)
  expect_equal(unname(out$cov), unname(hess$cov))
  # And 8 reflections place the normal where the mean is: the mode, here.
  expect_equal(out$details$normalcheck$centredraws, 8L)
  expect_equal(out$details$normalcheck$centre, rep(0, 6), tolerance = 1e-10)
  drawn <- ctsem:::.ctOptimDrawSamples(out, 'imis', list(), rep(0, 6), 50,
    lpg = lpg)
  expect_identical(drawn$draws, 'normal')
  expect_equal(dim(drawn$samples), c(50L, 6L))
})

test_that("checked centres its normal on the mean, not the mode", {
  # A Gaussian of 60 parameters whose normal has the widths wrong by up to 60%:
  # with far fewer draws than parameters the slopes are unresolved, and the
  # centre still does not move -- a draw and its reflection give the intercept
  # exactly, whatever the slopes.
  set.seed(8)
  k <- 60
  m <- stats::rnorm(k)
  prec <- exp(stats::runif(k, -1, 1))
  gb <- function(X) list(value = -0.5 * colSums(prec * (t(X) - m)^2),
    gradient = -t(prec * (t(X) - m)))
  check <- ctsem:::.ctOptimNormalCheck(gb, m,
    diag(exp(stats::runif(k, log(0.6), log(1.6))) / prec))
  expect_equal(check$centre, m, tolerance = 1e-10)

  # Six log-gamma coordinates, p(x) proportional to exp(a x - e^x): skewed, mode
  # log(a), mean digamma(a), 0.18 posterior sds apart at a = 8 -- too little
  # skew to sample, enough to move a mean. The draws are placed nearer it.
  a <- 8
  lpg <- function(x) {
    out <- sum(a * x - exp(x))
    attr(out, 'gradient') <- a - exp(x)
    out
  }
  target <- list(lpg = lpg, mode = rep(log(a), 6))
  set.seed(4)
  run <- .checkedRun(target, 'checked')
  expect_identical(run$draws, 'normal')
  sd <- sqrt(trigamma(a))
  frommode <- abs(log(a) - digamma(a)) / sd
  fromcentre <- mean(abs(run$check$centre - digamma(a)) / sd)
  expect_lt(fromcentre, 0.8 * frommode)
  expect_equal(colMeans(run$samples), run$check$centre, tolerance = 0.1)
})

test_that("checked's centre does not read a curvature it leaves unnamed as a shift", {
  # A five-dimensional shell of radius 0.05 under a N((0.5, 0, ...), I) prior,
  # mixed with 25 standard normal coordinates by a scaled rotation: too slight
  # in any one parameter's variance for the check to name, and the even part
  # of its residuals is the shell's curvature. Read as a shift it moved the
  # centre the full three sds, where the mode was within 0.02 posterior sds of
  # the mean; taken again about the moved centre it does not fall, and the
  # step goes.
  set.seed(11)
  mix <- qr.Q(qr(matrix(rnorm(900), 30))) %*% diag(exp(rnorm(30, 0, 0.5)))
  inv <- solve(mix)
  r0 <- 0.05; tau <- 0.01
  m0 <- c(0.5, 0, 0, 0, 0)
  lpg <- function(theta) {
    z <- as.numeric(inv %*% theta)
    x <- z[1:5]
    y <- z[6:30]
    r <- sqrt(sum(x^2))
    out <- -(r - r0)^2 / (2 * tau^2) - sum((x - m0)^2) / 2 - sum(y^2) / 2
    attr(out, 'gradient') <- as.numeric(t(inv) %*% c(-(r - r0) / tau^2 * x / r - (x - m0), -y))
    out
  }
  mode <- optim(as.numeric(mix %*% c(r0, rep(0, 29))), function(x) -lpg(x),
    function(x) -attr(lpg(x), 'gradient'), method = 'BFGS',
    control = list(reltol = 1e-15, maxit = 5000))$par
  hess <- suppressWarnings(ctsem:::ctOptimComputeUncertainty(mode,
    list(nsubjects = 1), sm = NULL, lpgFunc = lpg, uncertainty = 'hessian',
    finishsamples = 20, verbose = 0))
  gb <- attr(ctsem:::.ctImisPointGradbatch(lpg), 'gradbatch')
  sds <- sqrt(diag(hess$cov))
  moved <- vapply(1:5, function(s) {
    set.seed(s)
    check <- ctsem:::.ctOptimNormalCheck(gb, mode, hess$cov)
    max(abs(check$centre - mode) / sds)
  }, numeric(1))
  expect_lt(max(moved), 0.2)
})

test_that("checked samples a ring the Hessian's normal makes six times too wide", {
  # A shell of radius 0.05 under a N((0.5, 0), I) prior, rotated 45 degrees:
  # the direction around it is the prior's, the Hessian's curvature there a
  # tenth of a unit, and the normal's sd 0.24 against the ring's 0.037.
  r0 <- 0.05; tau <- 0.01
  rot <- matrix(c(1, 1, -1, 1), 2) / sqrt(2)
  lp <- function(X) {
    uv <- X %*% rot
    r <- sqrt(rowSums(uv^2))
    -(r - r0)^2 / (2 * tau^2) - ((uv[, 1] - 0.5)^2 + uv[, 2]^2) / 2
  }
  grad <- function(x) {
    uv <- as.numeric(t(rot) %*% x)
    r <- sqrt(sum(uv^2))
    as.numeric(rot %*% (-(r - r0) / tau^2 * uv / r - c(uv[1] - 0.5, uv[2])))
  }
  mode <- optim(as.numeric(rot %*% c(r0, 0)), function(x) -lp(matrix(x, 1)),
    function(x) -grad(x), method = 'BFGS', control = list(reltol = 1e-15))$par
  g <- seq(-0.12, 0.12, length.out = 401)
  target <- .checkedTarget(lp, grad, mode, as.matrix(expand.grid(g, g)))
  set.seed(1)
  hess <- .checkedRun(target, 'hessian')
  expect_gt(min(sqrt(diag(hess$uncertaintyfit$cov)) / target$sd), 5)
  set.seed(1)
  checked <- .checkedRun(target, 'checked')
  expect_equal(ncol(checked$check$directions), 2L)
  expect_identical(checked$draws, 'imis')
  expect_equal(checked$uncertaintyfit$details$importance_sampling$directions, 2L)
  expect_equal(sqrt(diag(checked$uncertaintyfit$cov)) / target$sd, c(1, 1),
    tolerance = 0.15)
})

test_that("checked's centre leaves the coordinates beside a shell where they are", {
  # A four-dimensional shell of radius 0.05 under a N((0.5, 0, 0, 0), I) prior,
  # with six exactly normal coordinates beside it: the check samples the shell,
  # and the centre's step in the rest must not borrow the shell's intercept.
  # Solved in every direction and projected afterwards it moved the six by up
  # to 10 sds, whose mean is 0.
  r0 <- 0.05; tau <- 0.01
  lpx <- function(x) {
    r <- sqrt(sum(x^2))
    -(r - r0)^2 / (2 * tau^2) - sum((x - c(0.5, 0, 0, 0))^2) / 2
  }
  gx <- function(x) {
    r <- sqrt(sum(x^2))
    -(r - r0) / tau^2 * x / r - (x - c(0.5, 0, 0, 0))
  }
  lpg <- function(th) {
    out <- lpx(th[1:4]) - sum(th[5:10]^2) / 2
    attr(out, 'gradient') <- c(gx(th[1:4]), -th[5:10])
    out
  }
  gb <- function(X) list(value = apply(X, 1, function(x) as.numeric(lpg(x))),
    gradient = t(apply(X, 1, function(x) attr(lpg(x), 'gradient'))))
  mode <- c(optim(c(r0, 0, 0, 0), function(x) -lpx(x), function(x) -gx(x),
    method = 'BFGS', control = list(reltol = 1e-15))$par, rep(0, 6))
  hess <- suppressWarnings(ctsem:::ctOptimComputeUncertainty(mode,
    list(nsubjects = 1), sm = NULL, lpgFunc = lpg, uncertainty = 'hessian',
    finishsamples = 20, verbose = 0))
  set.seed(11)
  check <- ctsem:::.ctOptimNormalCheck(gb, mode, hess$cov)
  expect_gte(ncol(check$directions), 1L)
  expect_equal(check$centre[5:10], rep(0, 6), tolerance = 1e-8)
})

test_that("checked samples a direction the Hessian gives no curvature", {
  # A ring of radius 0.05 under a N(0, I) prior centred where the ring is: at
  # its radial mode the direction around it has exactly no curvature, the
  # covariance projects it out, and the Hessian's normal has no spread along it
  # at all. The routes that draw give it the width the log density falls half
  # a unit over, and 'checked' then samples the ring.
  r0 <- 0.05; tau <- 0.01
  lp <- function(X) {
    r <- sqrt(rowSums(X^2))
    -(r - r0)^2 / (2 * tau^2) - rowSums(X^2) / 2
  }
  grad <- function(x) {
    r <- sqrt(sum(x^2))
    -(r - r0) / tau^2 * x / r - x
  }
  rs <- r0 / (1 + tau^2)
  mode <- c(rs, 0)
  hessian <- diag(c(-1 / tau^2 - 1, (r0 - rs) / (tau^2 * rs) - 1))
  g <- seq(-0.12, 0.12, length.out = 401)
  target <- .checkedTarget(lp, grad, mode, as.matrix(expand.grid(g, g)))
  set.seed(1)
  hess <- suppressWarnings(.checkedRun(target, 'hessian', hessian = hessian))
  expect_lt(sqrt(hess$uncertaintyfit$cov[2, 2]) / target$sd[2], 0.01)
  set.seed(1)
  checked <- suppressWarnings(.checkedRun(target, 'checked', hessian = hessian))
  details <- checked$uncertaintyfit$details
  widened <- (details$proposal_covariance %||% details$covariance)$widened
  expect_equal(widened$n, 1L)
  expect_equal(abs(widened$directions[, 1]), c(0, 1), tolerance = 1e-8)
  expect_gte(ncol(checked$check$directions), 1L)
  expect_equal(sqrt(diag(checked$uncertaintyfit$cov)) / target$sd, c(1, 1),
    tolerance = 0.15)
  # The projection dropped the tangent coordinate, and nothing is held once it
  # has a width: the intervals are read against what is held, so the summary
  # does not blank a coordinate the draws spread along.
  covariance <- details$proposal_covariance %||% details$covariance
  expect_equal(covariance$nullMass, c(0, 1), tolerance = 1e-8)
  expect_equal(covariance$heldMass, c(0, 0), tolerance = 1e-8)
  expect_equal(ctsem:::.ctBackendHeldMass(details), c(0, 0), tolerance = 1e-8)
  expect_equal(ctsem:::.ctBackendHeldMass(hess$uncertaintyfit$details), c(0, 1),
    tolerance = 1e-8)
})

test_that("checked corrects the funnel a variance from six observations makes", {
  # y ~ N(mu, sigma) from six values, theta = (mu, log sigma). The mode sits
  # where sigma is smallest, so the normal is narrow in mu and its mean of
  # log sigma low by half a posterior sd.
  y <- c(-0.3, 1.2, 0.4, 2.1, 0.9, 1.5)
  lp <- function(X) {
    s <- exp(X[, 2])
    -length(y) * X[, 2] - (colSums(outer(y, X[, 1], '-')^2)) / (2 * s^2) -
      X[, 1]^2 / 200 - X[, 2]^2 / 8
  }
  grad <- function(x) {
    s2 <- exp(2 * x[2])
    c(sum(y - x[1]) / s2 - x[1] / 100,
      -length(y) + sum((y - x[1])^2) / s2 - x[2] / 4)
  }
  mode <- optim(c(1, 0), function(x) -lp(matrix(x, 1)), function(x) -grad(x),
    method = 'BFGS', control = list(reltol = 1e-14))$par
  target <- .checkedTarget(lp, grad, mode, as.matrix(expand.grid(
    seq(-2, 4, length.out = 401), seq(-1.5, 1.8, length.out = 401))))
  set.seed(2)
  hess <- .checkedRun(target, 'hessian')
  expect_gt(max(abs(colMeans(hess$samples) - target$mean) / target$sd), 0.4)
  set.seed(2)
  checked <- .checkedRun(target, 'checked')
  expect_gte(ncol(checked$check$directions), 1L)
  expect_lt(max(abs(colMeans(checked$samples) - target$mean) / target$sd), 0.15)
  expect_equal(sqrt(diag(checked$uncertaintyfit$cov)) / target$sd, c(1, 1),
    tolerance = 0.15)
})

test_that("checked keeps the normal, and says so, where sampling cannot correct it", {
  # A four-dimensional shell of radius 0.05: the Hessian's normal is many
  # times too wide around it, and one round of 300 draws cannot place it.
  d <- 4; r0 <- 0.05; tau <- 0.01; e1 <- c(1, rep(0, d - 1))
  lpg <- function(x) {
    r <- sqrt(sum(x^2))
    out <- -(r - r0)^2 / (2 * tau^2) - sum((x - 0.5 * e1)^2) / 2
    attr(out, 'gradient') <- -(r - r0) / tau^2 * x / r - (x - 0.5 * e1)
    out
  }
  mode <- optim(r0 * e1, function(x) -lpg(x), function(x) -attr(lpg(x), 'gradient'),
    method = 'BFGS', control = list(reltol = 1e-15))$par
  set.seed(5)
  out <- suppressMessages(ctsem:::ctOptimComputeUncertainty(mode, list(nsubjects = 1),
    sm = NULL, lpgFunc = lpg, uncertainty = 'checked', finishsamples = 20, verbose = 0))
  expect_gt(ncol(out$details$normalcheck$directions), 0L)
  expect_true(all(2:d %in% out$details$normalcheck$parameters))
  expect_warning(drawn <- ctsem:::.ctOptimDrawSamples(out, 'imis',
    list(imisBudget = 300), mode, 200, lpg = ctsem:::.ctImisPointGradbatch(lpg),
    parnames = paste0('p', seq_len(d))), "unreliable for p")
  expect_identical(drawn$draws, 'normal')
  expect_equal(drawn$uncertaintyfit$cov, out$cov)
  expect_false(drawn$uncertaintyfit$details$importance_sampling$used)
})

test_that("the sampler's metric is the Hessian rescaled to the widths the check measured", {
  set.seed(6)
  P <- crossprod(matrix(rnorm(25), 5)) + diag(5)
  H <- -P
  eig <- eigen(solve(P), symmetric = TRUE)
  A <- eig$vectors %*% diag(sqrt(eig$values))
  v <- qr.Q(qr(matrix(rnorm(10), 5)))
  widths <- c(0.01, 4)
  fit <- list(uncertainty = list(hessian = H, details = list(normalcheck =
    list(directions = A %*% v, widths = widths))))
  rescaled <- -ctsem:::.ctBackendSampleMetric(fit)
  # Where the Hessian's normal is N(0, I), the precision is the identity but
  # 1 / width along each direction the check named.
  expect_equal(t(A) %*% rescaled %*% A, diag(5) + v %*% diag(1 / widths - 1) %*% t(v))
  expect_identical(rescaled, t(rescaled))
  # Nothing named, nothing rescaled: the engine reads the fit's own Hessian.
  fit$uncertainty$details$normalcheck$widths <- numeric()
  expect_null(ctsem:::.ctBackendSampleMetric(fit))
})

test_that("uncertainty = 'surrogate' is refused by name", {
  expect_error(ctFitUncertainty(list(), 'surrogate'), "removed in ctsem 4.0.0")
  expect_error(ctOptimUncertainty(list(), 'surrogate'), "uncertainty = 'checked'")
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
  expect_identical(eval(formals(ctFitUncertainty)$uncertainty)[1], 'checked')
  expect_false('surrogate' %in% eval(formals(ctFitUncertainty)$uncertainty))
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

# Which flat directions are flat in the *likelihood*, not just in a
# numerically differentiated approximation to it.
#
# The curvature along a direction the data does not determine is not a property
# of the model and the data. Measured on a two-latent model fitted to noise,
# walking one diffusion correlation out along its flat ray: the log likelihood
# is -207.01897 at every raw value from -6 to -20, while the smallest relative
# eigenvalue falls from 1.6e-08 to 7.1e-16 and then turns negative from
# rounding. So which side of `.ctFlatDirectionRtol()` the direction landed on
# was decided by where the optimiser stopped, and the same fit with its
# predicted-gain rule on and off got opposite diagnoses.
#
# The screen asks the likelihood instead, against the likelihood-ratio bound
# (Raue et al. 2009), which is a statistical quantity rather than a tolerance
# on an approximation. Tested here on functions whose answer is known, because
# the whole question is which functions it calls flat.
test_that("a flat direction is confirmed against the likelihood, not the curvature", {
  screen <- ctsem:::.ctOptimFlatDirectionScreen
  at <- c(0, 0)
  bar <- stats::qchisq(0.95, 1) / 2

  # Genuinely flat along the second coordinate, and the curvature agrees.
  flat <- screen(diag(c(1, 1e-14)), function(x) -0.5 * x[1]^2, at)
  expect_equal(flat$flat, c(FALSE, TRUE))
  expect_equal(flat$change[2], 0)
  expect_equal(flat$bar, bar)

  # A curvature that says flat where the likelihood is not. This is the case
  # the screen exists to refuse, and refusing it is what makes the screen safe
  # to act on: it can only ever remove a direction from the identified
  # subspace, so a false positive here would invent a missing interval.
  lying <- screen(diag(c(1, 1e-14)), function(x) -0.5 * sum(x^2), at)
  expect_equal(lying$flat, c(FALSE, FALSE))
  expect_gt(lying$change[2], bar)

  # A direction the likelihood *rises* along is not flat -- it is one the
  # optimiser has not finished with. It arrives as a candidate because a
  # rounding-negative eigenvalue is below any threshold, and it must not be
  # confirmed, which is why the change is measured as a magnitude rather than
  # as a drop.
  rising <- screen(diag(c(1, -1e-16)),
    function(x) -0.5 * x[1]^2 + 3 * abs(x[2]), at)
  expect_equal(rising$flat, c(FALSE, FALSE))
  expect_gt(rising$change[2], bar)

  # And the usual fit, which has no flat direction: nothing to ask about, so
  # nothing is evaluated. This is what keeps the screen free on the fits that
  # are the overwhelming majority -- it returns before the first call.
  calls <- 0L
  counted <- function(x) { calls <<- calls + 1L; -0.5 * sum(x^2) }
  expect_null(screen(diag(c(1, 0.5)), counted, at))
  expect_equal(calls, 0L)

  # A point the model cannot evaluate is not evidence of flatness either.
  broken <- screen(diag(c(1, 1e-14)),
    function(x) if (abs(x[2]) > 1e-8) NaN else -0.5 * x[1]^2, at)
  expect_equal(broken$flat, c(FALSE, FALSE))
  expect_true(is.na(broken$change[2]))
})

test_that("the confirmed directions come out of the covariance, and only those", {
  eig <- eigen(diag(c(1, 0.5)), symmetric = TRUE)
  # The mask only ever removes. With none set this is the eigenvalue rule it
  # has always been, which is what makes it safe to leave on everywhere.
  plain <- ctsem:::.ctOptimIdentifiedInverse(diag(c(1, 0.5)))
  expect_equal(plain$nnull, 0L)
  masked <- ctsem:::.ctOptimIdentifiedInverse(diag(c(1, 0.5)), eig = eig,
    flat = c(FALSE, TRUE))
  expect_equal(masked$nnull, 1L)
  expect_equal(masked$nullMass, c(0, 1))

  # And `ctOptimCovFromHessian()` takes the projection branch on measured
  # evidence even where the eigenvalue alone would have let `solve()` through.
  # That ordering matters for the reason the existing comment there gives:
  # whether `solve()` succeeds on a nearly singular matrix is settled by
  # rounding, so a repair reached only on failure is reached only sometimes.
  info <- diag(c(1, 1e-10))
  loose <- suppressWarnings(suppressMessages(
    ctsem:::ctOptimCovFromHessian(-info, warn = FALSE)))
  expect_equal(attr(loose, 'ctOptimCovFromHessian')$method, 'solve')
  confirmed <- ctsem:::.ctOptimFlatDirectionScreen(info,
    function(x) -0.5 * x[1]^2, c(0, 0))
  expect_true(any(confirmed$flat))
  screened <- suppressWarnings(suppressMessages(
    ctsem:::ctOptimCovFromHessian(-info, warn = FALSE, screen = confirmed)))
  diagnostics <- attr(screened, 'ctOptimCovFromHessian')
  expect_equal(diagnostics$method, 'nullprojection')
  expect_equal(diagnostics$profileFlatDirections, 1L)
  expect_lt(max(diagnostics$profileChange), diagnostics$profileBar)
})

# A ridge that is flat but curved in raw coordinates, which is what the
# stan-julia parity fixture's ten correlations lie on. A straight slice leaves a
# curved ridge, and how fast depends on where along it the optimiser stopped --
# measured there at 17.8 nats at four raw units from one stopping point and 0.68
# from another, 4.4e-04 nats higher on the same ridge. So a rung that drops past
# the bar is followed back to the ridge before it is judged, and one flat side is
# enough.
#
# f(x, y) = -(y - x^2 / 2)^2 is exactly flat along the parabola y = x^2 / 2 and
# at the origin its curvature along x is exactly zero, so x is the candidate.
# The straight walk to x = 4 drops 64 nats; one Newton step in y finds the
# ridge again.
.parabolic_ridge <- function(p) {
  value <- -(p[2] - p[1]^2 / 2)^2
  attr(value, 'gradient') <- c(2 * (p[2] - p[1]^2 / 2) * p[1],
    -2 * (p[2] - p[1]^2 / 2))
  value
}

test_that("a curved ridge is followed, and confirmed flat along it", {
  info <- diag(c(0, 2))
  bar <- stats::qchisq(0.95, 1) / 2
  expect_gt(-as.numeric(.parabolic_ridge(c(4, 0))), bar)
  found <- ctsem:::.ctOptimFlatDirectionScreen(info, .parabolic_ridge, c(0, 0))
  flat <- which(found$flat)
  expect_length(flat, 1L)
  # The flat direction is the x axis.
  expect_equal(abs(found$eig$vectors[, flat]), c(1, 0))
  expect_lt(found$change[flat], bar)
  # And a correction needs the gradient: the same function without one is a
  # straight slice, which this ridge refuses.
  bare <- ctsem:::.ctOptimFlatDirectionScreen(info,
    function(p) as.numeric(.parabolic_ridge(p)), c(0, 0))
  expect_false(any(bare$flat))
  expect_gt(max(bare$change, na.rm = TRUE), bar)
})

test_that("one flat side settles it; a side that cannot be evaluated says nothing", {
  bar <- stats::qchisq(0.95, 1) / 2
  info <- diag(c(1, 1e-14))
  # Flat for positive y, falling steeply for negative: a curve of constant
  # likelihood runs four units out on one side, which is what non-identification
  # needs. Which side the eigenvector points to is arbitrary, so the likelihood
  # is flat on whichever the other one is not.
  halfflat <- function(p) -0.5 * p[1]^2 - 10 * min(p[2], 0)^2
  found <- ctsem:::.ctOptimFlatDirectionScreen(info, halfflat, c(0, 0))
  expect_equal(found$flat, c(FALSE, TRUE))
  expect_equal(found$change[2], 0)
  # Steep on both sides is refused, whichever side comes first.
  steep <- function(p) -0.5 * p[1]^2 - 10 * p[2]^2
  expect_false(any(ctsem:::.ctOptimFlatDirectionScreen(info, steep,
    c(0, 0))$flat))
  # A side the model cannot evaluate is not evidence; the other side still is.
  onesided <- function(p) if (p[2] < -1e-8) NaN else -0.5 * p[1]^2
  found <- ctsem:::.ctOptimFlatDirectionScreen(info, onesided, c(0, 0))
  expect_equal(found$flat, c(FALSE, TRUE))
})

test_that("candidates reach curvature a stopping point has not yet decayed", {
  # The parity fixture's ridge sat at 1.7e-08 of the sharpest curvature at its
  # earliest stopping point, above the eigenvalue rule's 1e-8, while its first
  # identified direction held at 1.2e-06 everywhere. The gate is between them.
  flatline <- function(p) -0.5 * p[1]^2
  expect_true(ctsem:::.ctOptimFlatDirectionScreen(diag(c(1, 3e-8)),
    flatline, c(0, 0))$flat[2])
  # And what is above it costs nothing: nothing is evaluated.
  calls <- 0L
  counted <- function(p) { calls <<- calls + 1L; -0.5 * p[1]^2 }
  expect_null(ctsem:::.ctOptimFlatDirectionScreen(diag(c(1, 3e-7)), counted,
    c(0, 0)))
  expect_equal(calls, 0L)
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

test_that("the projection warning says it once, and fits in a warning", {
  info <- .mixed_null_information()
  text <- tryCatch(ctsem:::ctOptimCovFromHessian(-info),
    warning = function(w) conditionMessage(w))

  # It used to say the same thing three times -- solve was skipped, directions
  # were projected out, they have no reported spread -- and the total ran past
  # R's 1000-byte warning cap, so the end was cut off mid-word.
  expect_lt(nchar(text, type = "bytes"), 1000L)
  expect_false(grepl("not attempted", text, fixed = TRUE))
  expect_false(grepl("projected out before inversion", text, fixed = TRUE))
  # Said once, with both counts: directions dropped and parameters affected.
  expect_match(text, "1 direction\\(s\\) with no curvature left out")
  expect_match(text, "2 parameter\\(s\\) along them")

  # The audit trail is unchanged -- the steps are still all on the object, it
  # is only the warning that is a phrase.
  cov <- suppressWarnings(suppressMessages(
    ctsem:::ctOptimCovFromHessian(-info, warn = FALSE)))
  steps <- attr(cov, 'ctOptimCovFromHessian')$repairSteps
  expect_true(any(grepl("not attempted", steps, fixed = TRUE)))
  expect_true(any(grepl("projected out before inversion", steps, fixed = TRUE)))
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
  # `.ctParticleCorrect()` and `.ctLaplaceCorrect()` rely on when they pre-scale
  # their own proposal and pass 1: the weighted covariance recovers the target.
  unscaled <- run(Sigma, 1)
  expect_equal(unscaled$covariance, Sigma, tolerance = 0.1)
})

# IS-importance-sampling-2026-09-06.md: an importance sampler asked to sample a
# raw direction with no curvature at all is asking it to sample a density that
# is not one -- the weights have infinite variance and the effective sample
# size never converges (measured there: 51,000 evaluations, ESS oscillating
# between 1.1 and 64.8 against a target of 100, and standard errors of 873 and
# 1214 on parameters the Hessian puts at 1.5 and 0.8). `.ctImisSubspace()` and
# `.ctImisRun()` are the repair: run IMIS in the whitened eigen-coordinates of
# the directions a covariance actually has curvature in, and hold the rest at
# the estimate rather than sampling them at all.
#
# These tests are deterministic and need no fit: the target below is exactly
# Gaussian in three dimensions and exactly flat in a fourth, which is what a
# raw correlation with no identified individual differences behind it looks
# like, and the closed form is the reference rather than a stored number.
.imis_subspace_fixture <- function() {
  Q <- qr.Q(qr(matrix(c(
    0.5, -0.3, 0.7, 0.2, 0.6, 0.5, -0.2, -0.4,
    -0.3, 0.7, 0.1, 0.6, 0.4, -0.4, 0.8, 0.1), 4, 4)))
  eigvals <- c(4, 1.5, 0.6, 0) # the fourth direction: no curvature at all
  info <- Q %*% diag(1 / ifelse(eigvals == 0, Inf, eigvals)) %*% t(Q)
  info[!is.finite(info)] <- 0
  # The covariance `.ctOptimIdentifiedInverse()` would hand back for this
  # information matrix: zero variance along the null direction, the true
  # covariance along the other three -- not a ridge-floored stand-in.
  cov <- Q %*% diag(eigvals) %*% t(Q)
  centre <- c(1.2, -0.5, 0.3, 2.0)
  target <- function(x) -0.5 * as.numeric(t(x - centre) %*% info %*% (x - centre))
  list(Q = Q, cov = cov, centre = centre, target = target)
}

test_that("a covariance with no flat direction is passed to imis_is unchanged", {
  skip_if_not_installed('mvtnorm'); skip_if_not_installed('diagis')
  skip_if_not_installed('gridExtra'); skip_if_not_installed('ggplot2')

  fullrank <- diag(c(1, 2, 3, 4))
  expect_null(ctsem:::.ctImisSubspace(fullrank))

  target <- function(x) -0.5 * sum(x^2 / c(1, 2, 3, 4))
  direct <- withr::with_seed(5, ctsem:::imis_is(target, mu_hat = rep(0, 4),
    Sigma_hat = fullrank, max_iter = 2, scale_init = 1.5, tail_scale = 1.2,
    df = Inf, target_ess = 1e9, n_batch = 200, cl = NA, finishsamples = 30,
    verbose = FALSE, diag_plots = FALSE))
  viarun <- withr::with_seed(5, ctsem:::.ctImisRun(target, centre = rep(0, 4),
    cov = fullrank, max_iter = 2, scale_init = 1.5, tail_scale = 1.2,
    df = Inf, target_ess = 1e9, n_batch = 200, cl = NA, finishsamples = 30,
    verbose = FALSE, diag_plots = FALSE))
  expect_identical(direct$full_theta, viarun$full_theta)
  expect_null(attr(viarun, 'subspace'))
})

test_that("IMIS in the identified subspace holds the null direction at the estimate", {
  skip_if_not_installed('mvtnorm'); skip_if_not_installed('diagis')
  skip_if_not_installed('gridExtra'); skip_if_not_installed('ggplot2')

  fx <- .imis_subspace_fixture()
  subspace <- ctsem:::.ctImisSubspace(fx$cov)
  expect_equal(subspace$k, 3L)
  expect_equal(subspace$nnull, 1L)

  result <- withr::with_seed(11, ctsem:::.ctImisRun(fx$target, centre = fx$centre,
    cov = fx$cov, max_iter = 50, scale_init = 1.5, tail_scale = 1.2, df = Inf,
    target_ess = 100, n_batch = 1000, cl = NA, finishsamples = 300,
    verbose = FALSE, diag_plots = FALSE))
  expect_equal(attr(result, 'subspace')$nnull, 1L)
  expect_true(result$ess >= 100)
  # Reached the target from far fewer evaluations than the unwhitened run
  # needs on a genuinely flat direction (13x to 26x fewer, measured in the IS
  # note); one iteration's worth (1000) here is generous rather than a tight
  # bound, since the point is convergence, not the exact count.
  expect_lte(nrow(result$full_theta), 2000L)

  nulldir <- fx$Q[, 4]
  # No draw ever moves off the affine subspace through `centre` that the three
  # kept eigenvectors span -- not "small", exactly zero up to the eigenbasis's
  # own floating-point orthogonality (~1e-15 at this scale).
  offset <- sweep(result$full_theta, 2, fx$centre, '-')
  expect_lt(max(abs(offset %*% nulldir)), 1e-8)
  expect_lt(abs(as.numeric(t(nulldir) %*% result$covariance %*% nulldir)), 1e-8)
  expect_equal(as.numeric(t(result$mean) %*% nulldir),
    as.numeric(t(fx$centre) %*% nulldir), tolerance = 1e-8)

  # And the three identified directions recover the true covariance -- the
  # closed-form reference, not a second sampler's opinion of it.
  Vid <- fx$Q[, 1:3]
  expect_equal(t(Vid) %*% result$covariance %*% Vid, t(Vid) %*% fx$cov %*% Vid,
    tolerance = 0.2)
})

test_that("whitening carries a density's 'batch' attribute through the map", {
  skip_if_not_installed('mvtnorm'); skip_if_not_installed('diagis')
  skip_if_not_installed('gridExtra'); skip_if_not_installed('ggplot2')

  fx <- .imis_subspace_fixture()
  calls_scalar <- 0L
  calls_batch <- 0L
  target <- function(x) { calls_scalar <<- calls_scalar + 1L; fx$target(x) }
  attr(target, 'batch') <- function(X) {
    calls_batch <<- calls_batch + 1L
    apply(X, 1, fx$target)
  }
  result <- withr::with_seed(11, ctsem:::.ctImisRun(target, centre = fx$centre,
    cov = fx$cov, max_iter = 5, scale_init = 1.5, tail_scale = 1.2, df = Inf,
    # Unreachable, so every one of the six iterations (0:5) runs and the count
    # is exact rather than "at least one".
    target_ess = 1e9, n_batch = 50, cl = NA, finishsamples = 20,
    verbose = FALSE, diag_plots = FALSE))
  expect_identical(calls_scalar, 0L)
  expect_identical(calls_batch, 6L)
  expect_true(is.finite(result$ess))
})

test_that("checked's sampler can weigh the sampled directions by their marginal", {
  skip_if_not_installed('mvtnorm'); skip_if_not_installed('diagis')
  skip_if_not_installed('gridExtra'); skip_if_not_installed('ggplot2')
  # s ~ N(0, 1) and eight coordinates c | s ~ N(0, e^(a s)). The mode is at
  # s = -4a: held at c = 0 there, s has the conditional N(-4a, 1); its
  # marginal is N(0, 1), the volume the c's spread adds as s grows.
  a <- 0.1
  k <- 8L
  lp <- function(x) -x[1]^2 / 2 - sum(x[-1]^2) * exp(-a * x[1]) / 2 - k * a * x[1] / 2
  attr(lp, 'batch') <- function(X) apply(as.matrix(X), 1, lp)
  mode <- c(-k * a / 2, rep(0, k))
  cov <- diag(c(1, rep(exp(-k * a^2 / 2), k)))
  run <- function(complement) withr::with_seed(5, ctsem:::.ctImisRun(lp, centre = mode,
    cov = cov, directions = matrix(c(1, rep(0, k))), widths = 1, complement = complement,
    max_iter = 50, scale_init = 1.5, tail_scale = 1.2, df = 5, target_ess = 1000,
    n_batch = 500, cl = NA, finishsamples = 1000, paths = FALSE, verbose = FALSE,
    diag_plots = FALSE))
  held <- run(FALSE)
  integrated <- run(TRUE)
  expect_lt(abs(held$mean[1] + k * a / 2), 0.12)
  expect_lt(abs(integrated$mean[1]), 0.12)
  # The others' answer is the normal's either way.
  expect_equal(diag(integrated$covariance)[-1], diag(cov)[-1])
})

test_that(".ctOptimImisDraws() reports which directions the subspace held", {
  skip_if_not_installed('mvtnorm'); skip_if_not_installed('diagis')
  skip_if_not_installed('gridExtra'); skip_if_not_installed('ggplot2')

  fx <- .imis_subspace_fixture()
  drawn <- withr::with_seed(11, ctsem:::.ctOptimImisDraws(fx$target,
    centre = fx$centre, cov = fx$cov, finishsamples = 100,
    remedy = "test remedy", nbatch = 1000, target_ess = 100, maxiter = 50,
    scaleInit = 1.5, tailScale = 1.2, diagPlots = FALSE))
  expect_equal(drawn$subspace$nnull, 1L)
  expect_true(drawn$ess >= 100)
  expect_equal(ncol(drawn$samples), 4L)
})

test_that("importance sampling records the Pareto k of its weights and warns when it is high", {
  skip_if_not_installed('mvtnorm'); skip_if_not_installed('diagis')
  skip_if_not_installed('gridExtra'); skip_if_not_installed('ggplot2')
  skip_if_not_installed('loo')

  # Two closed-form targets for one standard normal proposal (`imisDf = Inf`;
  # the default t would bound these weights). A normal three times wider than
  # the proposal gives weights p/q proportional to exp(|x|^2 (1 - 1/9) / 2),
  # whose tail under the proposal is an exact power law with k = 1 - 1/9 =
  # 0.89: no variance, however many draws. A normal narrower than the proposal
  # gives bounded weights and k well below 0.5. One batch of 2000 each, so
  # nothing adapts; the densities carry no gradient, so no path is searched.
  heavy <- function(x) sum(stats::dnorm(x, 0, 3, log = TRUE))
  light <- function(x) sum(stats::dnorm(x, 0, 0.8, log = TRUE))
  draw <- function(lpg) withr::with_seed(3, ctsem:::.ctOptimDrawSamples(
    list(cov = diag(2), details = list()), draws = 'imis',
    control = list(imisMaxIter = 0, isitersize = 2000, isESS = 1, imisDf = Inf),
    est = c(0, 0), finishsamples = 50, lpg = lpg, scaleInit = 1, tailScale = 1,
    tailremedy = 'Sample it instead.'))

  expect_warning(bad <- draw(heavy), 'Pareto k .*Sample it instead')
  kbad <- bad$uncertaintyfit$details$importance_sampling$pareto_k
  expect_gt(kbad, 0.7)

  expect_no_warning(good <- draw(light))
  expect_lt(good$uncertaintyfit$details$importance_sampling$pareto_k, 0.5)

  # The shared route .ctLaplaceCorrect() and .ctParticleCorrect() take returns
  # it too, and warns with the caller's own remedy.
  expect_warning(drawn <- withr::with_seed(3, ctsem:::.ctOptimImisDraws(heavy,
    centre = c(0, 0), cov = diag(2), finishsamples = 50, remedy = 'Caller remedy.',
    nbatch = 2000, target_ess = 1, maxiter = 0, scaleInit = 1, tailScale = 1,
    diagPlots = FALSE)), 'Pareto k .*Caller remedy')
  expect_equal(drawn$k, kbad)
})

# A variance-like tail along a curved ridge, with exact moments: x1 = -1.1 -
# exp(u), u ~ N(log .15, 1), and x2 = -0.7 (u - log .15) + sqrt(.51) e moving
# with it, so x2 ~ N(0, 1) exactly; two more coordinates N(0, 1). x1's sd is
# about six standard errors of the curvature at the mode, and its 2.5%
# quantile about 18 out, along a path on which x2 moves too -- the shape the
# variances of the bench models had, where a proposal built on the curvature
# reported half the width. Both a value batch and a gradient batch, as the
# julia route gives.
.imis_curved_fixture <- function() {
  mu0 <- log(0.15); a <- 0.7; b <- sqrt(1 - a^2)
  lp1 <- function(x) {
    y <- -1.1 - x[1]
    if (y <= 0) return(-1e100)
    u <- log(y); z2 <- (x[2] + a * (u - mu0)) / b
    -0.5 * (u - mu0)^2 - u - 0.5 * z2^2 - 0.5 * sum(x[-(1:2)]^2)
  }
  gr1 <- function(x) {
    y <- -1.1 - x[1]
    if (y <= 0) return(rep(0, length(x)))
    u <- log(y); z2 <- (x[2] + a * (u - mu0)) / b
    c((-(u - mu0) - 1 - z2 * a / b) * (-1 / y), -z2 / b, -x[-(1:2)])
  }
  dens <- function(x) lp1(x)
  attr(dens, 'batch') <- function(X) apply(as.matrix(X), 1, lp1)
  attr(dens, 'gradbatch') <- function(X) {
    X <- as.matrix(X)
    list(value = apply(X, 1, lp1), gradient = t(apply(X, 1, gr1)))
  }
  mode <- stats::optim(c(-1.2, 0, 0, 0), function(x) -lp1(x), function(x) -gr1(x),
    method = 'BFGS', control = list(reltol = 1e-14, maxit = 1000))$par
  cov <- solve(stats::optimHess(mode, function(x) -lp1(x), function(x) -gr1(x)))
  list(dens = dens, mode = mode, cov = cov,
    sd = sqrt((exp(1) - 1) * 0.15^2 * exp(1)),
    q025 = -1.1 - 0.15 * exp(stats::qnorm(0.975)))
}

test_that("importance sampling follows a curved tail the curvature at the mode does not show", {
  skip_if_not_installed('mvtnorm'); skip_if_not_installed('diagis')
  skip_if_not_installed('gridExtra'); skip_if_not_installed('ggplot2')

  fx <- .imis_curved_fixture()
  # uncertainty = 'is' at its defaults, over four seeds, as sd ratio and 2.5%
  # quantile error in true sds. Measured: 0.82-1.01 and 0.01-0.54 with the
  # path search, medians 0.86 and 0.37; without it 0.43-0.80 and 0.07-1.81,
  # medians 0.67 and 0.64.
  run <- function(dens) sapply(1:4, function(s) {
    r <- withr::with_seed(s, suppressWarnings(ctsem:::.ctOptimDrawSamples(
      list(cov = fx$cov, details = list()), draws = 'imis', control = list(),
      est = fx$mode, finishsamples = 2000, lpg = dens)))
    x1 <- r$samples[, 1]
    c(sdratio = stats::sd(x1) / fx$sd,
      q025 = abs(stats::quantile(x1, 0.025, names = FALSE) - fx$q025) / fx$sd,
      followed = 1 %in% r$uncertaintyfit$details$importance_sampling$paths$followed)
  })
  out <- run(fx$dens)
  expect_true(all(out['followed', ] == 1))
  expect_gt(stats::median(out['sdratio', ]), 0.8)
  expect_lt(stats::median(out['q025', ]), 0.45)
})

test_that("a posterior that has not fallen off far out is reported", {
  skip_if_not_installed('mvtnorm'); skip_if_not_installed('diagis')
  skip_if_not_installed('gridExtra'); skip_if_not_installed('ggplot2')

  # A unit bump on a floor: curvature at the mode, then flat for ever along
  # x1 -- improper, as a random-effect sd with no prior can be.
  lp1 <- function(x) log(0.1 + exp(-0.5 * x[1]^2)) - 0.5 * x[2]^2
  gr1 <- function(x) c(-x[1] * exp(-0.5 * x[1]^2) / (0.1 + exp(-0.5 * x[1]^2)), -x[2])
  dens <- function(x) lp1(x)
  attr(dens, 'batch') <- function(X) apply(as.matrix(X), 1, lp1)
  attr(dens, 'gradbatch') <- function(X) {
    X <- as.matrix(X)
    list(value = apply(X, 1, lp1), gradient = t(apply(X, 1, gr1)))
  }
  expect_warning(r <- withr::with_seed(1, ctsem:::.ctOptimDrawSamples(
    list(cov = diag(c(1.1, 1)), details = list()), draws = 'imis',
    control = list(imisMaxIter = 2), est = c(0, 0), finishsamples = 100,
    lpg = dens, tailremedy = 'Sample it instead.')),
    'had not fallen off 96 standard errors .*raw parameters 1 \\(above\\), 1 \\(below\\).*Sample it instead')
  expect_setequal(r$uncertaintyfit$details$importance_sampling$paths$reachesLimit, c(1L, -1L))
})
