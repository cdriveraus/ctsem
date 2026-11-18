# Smooth TI predictor terms, s() in TIpredNames (R/ctTIpredSmooth.R).
#
# The checks against independent answers come first and fit nothing: the
# basis's penalty is the integrated squared second derivative with a line as
# its null space, and the Fellner-Schall update reaches the REML smoothing
# parameter of a Gaussian mixed model, which has a closed form. Then one fit
# with a curve in age, recovered and reported.

test_that("TIpredNames takes a formula of names and smooths", {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    LAMBDA = diag(1), CINT = matrix("cint"), MANIFESTMEANS = matrix(0),
    TIpredNames = ~ s(age, k = 6) + female))
  expect_equal(m$TIpredNames, c("age_lin", paste0("age_s", 1:4), "female"))
  expect_equal(m$TIsmooths$age$penalised, paste0("age_s", 1:4))
  # A smooth reaches the parameters that vary between persons by default; a
  # linear predictor every free one, as before.
  pars <- m$pars[!is.na(m$pars$param), ]
  varies <- pars$indvarying %in% TRUE
  expect_true(all(pars$age_s1_effect[varies] == "TRUE"))
  expect_true(all(pars$age_s1_effect[!varies] == "FALSE"))
  expect_true(all(pars$female_effect == "TRUE"))
  # A cell names the covariate to take the smooth.
  m2 <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    LAMBDA = diag(1), CINT = matrix("cint"), MANIFESTMEANS = matrix(0),
    DRIFT = matrix("dr||||age"), TIpredNames = ~ s(age, k = 5)))
  dr <- m2$pars$param %in% "dr"
  expect_true(all(unlist(m2$pars[dr, paste0(m2$TIpredNames, "_effect")]) == "TRUE"))
  expect_error(suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    LAMBDA = diag(1), DRIFT = matrix("dr||||age=0.5"),
    TIpredNames = ~ s(age, k = 5))), "estimated whole")
  expect_error(suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    LAMBDA = diag(1), DRIFT = matrix("dr||||age_s1"),
    TIpredNames = ~ s(age, k = 5))), "name the covariate")
  expect_error(ctsem:::.ctTIpredParseFormula(~ s(age, k = 3)), "at least 4")
  expect_error(ctsem:::.ctTIpredParseFormula(~ s(age) + age), "both")
  expect_error(ctsem:::.ctTIpredParseFormula(~ log(age)), "covariate names")
})

test_that("the basis is centred, its line unpenalised and the rest an identity penalty", {
  set.seed(3)
  x <- runif(200)
  basis <- ctsem:::.ctTIsmoothBuild(x, 8, "x")
  X <- ctsem:::.ctTIsmoothEval(basis, x)
  expect_equal(unname(colMeans(X)), rep(0, 7), tolerance = 1e-10)
  # The first column is the centred covariate at unit spread.
  expect_equal(as.numeric(X[, 1]), as.numeric(scale(x)) * sqrt(200 / 199),
    tolerance = 1e-8)
  # The integrated squared second derivative, by a fine Riemann sum: zero on
  # the line, the identity (up to one scale) on the rest.
  fine <- seq(0, 1, length.out = 20001)
  B2 <- ctsem:::.ctTIsmoothBSpline(basis$knots, fine, 2L) %*% basis$transform
  S <- crossprod(B2) * (fine[2] - fine[1])
  expect_equal(unname(S / S[2, 2]), diag(c(0, rep(1, 6))), tolerance = 1e-3)
  # Beyond the data the curve continues as a line.
  out <- ctsem:::.ctTIsmoothEval(basis, c(-1, -0.5, 0))
  expect_equal(unname(out[2, ] - out[1, ]), unname(out[3, ] - out[2, ]),
    tolerance = 1e-8)
})

test_that("the Fellner-Schall update reaches the REML smoothing parameter", {
  # y = X0 beta + Z b + e, e ~ N(0, 1): a Gaussian mixed model, where the
  # quadratic model is exact and REML has a closed form.
  set.seed(3)
  n <- 200
  x <- runif(n)
  X0 <- cbind(1, x)
  basis <- ctsem:::.ctTIsmoothBuild(x, 8, "x")
  Z <- ctsem:::.ctTIsmoothEval(basis, x)[, -1]
  y <- 1 + 2 * x + sin(6 * x) + rnorm(n)
  C <- cbind(X0, Z)
  m <- ncol(Z)
  pen <- 2 + seq_len(m)
  penalty <- data.frame(term = "s(x)", param = "p", m = m, lambda = 10)
  penalty$index <- list(pen)
  info <- crossprod(C)
  info[cbind(pen, pen)] <- info[cbind(pen, pen)] + 10
  theta <- as.numeric(solve(info, crossprod(C, y)))
  gradient <- as.numeric(crossprod(C, y - C %*% theta)) - c(0, 0, rep(10, m)) * theta
  gfs <- ctsem:::.ctTIsmoothLambda(penalty, -info, theta, gradient)
  reml <- function(logtau2) {
    V <- diag(n) + exp(logtau2) * tcrossprod(Z)
    Vi <- solve(V)
    XtVi <- crossprod(X0, Vi)
    beta <- solve(XtVi %*% X0, XtVi %*% y)
    r <- y - X0 %*% beta
    -0.5 * (determinant(V)$modulus + determinant(XtVi %*% X0)$modulus +
      crossprod(r, Vi %*% r))
  }
  best <- stats::optimize(function(l) -reml(l), c(-15, 10), tol = 1e-10)$minimum
  expect_equal(gfs$lambda, exp(-best), tolerance = 1e-4)
  # And the quadratic model's optimum is the penalised fit at that lambda.
  info2 <- crossprod(C)
  info2[cbind(pen, pen)] <- info2[cbind(pen, pen)] + gfs$lambda
  expect_equal(gfs$point, as.numeric(solve(info2, crossprod(C, y))),
    tolerance = 1e-8)
  # With a parameter the data say nothing about, the information is singular,
  # as a large fit's is along its flat directions: the update reads the
  # directions with curvature and reaches the same smoothing parameter, and the
  # point stays put along the one without. `solve()` stopped the affect
  # pilot's fit there after an 11-hour optimisation.
  flat <- matrix(0, ncol(C) + 1, ncol(C) + 1)
  flat[seq_len(ncol(C)), seq_len(ncol(C))] <- info
  again <- ctsem:::.ctTIsmoothLambda(penalty, -flat, c(theta, 0.3), c(gradient, 0))
  expect_equal(again$lambda, gfs$lambda, tolerance = 1e-8)
  expect_equal(again$edf, gfs$edf, tolerance = 1e-8)
  expect_equal(again$point, c(gfs$point, 0.3), tolerance = 1e-8)
})

test_that("smooths are refused where they cannot be estimated", {
  dat <- data.frame(id = rep(1:3, each = 3), time = rep(0:2, 3), Y1 = rnorm(9),
    age = rep(c(20, 30, 40), each = 3))
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    LAMBDA = diag(1), CINT = matrix("cint"), MANIFESTMEANS = matrix(0),
    TIpredNames = ~ s(age, k = 4)))
  expect_error(suppressMessages(ctFit(dat, m, backend = "stan", fit = FALSE)),
    "backend='julia'")
  expect_error(suppressMessages(ctFit(dat, m, backend = "julia",
    optimize = FALSE, fit = FALSE)), "optimize=TRUE")
  within <- dat
  within$age[2] <- 21
  expect_error(suppressMessages(ctFit(within, m, backend = "julia", fit = FALSE)),
    "varies within persons")
  few <- dat
  expect_error(ctsem:::.ctTIsmoothData(m, few), "distinct values")
  # A person without a value is held at the value the imputation model
  # predicts for them; here, with nothing else known, their study's mean shrunk
  # towards everyone's, and everyone's where the study has no values.
  set.seed(1)
  persons <- data.frame(id = 1:12, study = rep(1:3, each = 4),
    age = c(21, 25, 33, NA, 61, 67, 72, 78, NA, NA, NA, NA))
  long <- persons[rep(1:12, each = 2), ]
  long$time <- rep(0:1, 12)
  long$Y1 <- rnorm(24)
  long$age[3] <- NA
  nested <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    LAMBDA = diag(1), CINT = matrix("cint"), MANIFESTMEANS = matrix(0),
    id = c("id", "study"), TIpredNames = ~ s(age, k = 4)))
  out <- suppressWarnings(expect_message(ctsem:::.ctTIsmoothData(nested, long,
    includeOutcome = FALSE), "5 persons"))
  basis <- out$model$TIsmooths$age$basis
  at <- function(a) as.numeric(ctsem:::.ctTIsmoothEval(basis, a))
  row <- function(i) as.numeric(out$data[out$data$id == i,
    out$model$TIsmooths$age$columns][1L, ])
  # A value on any of a person's rows is theirs.
  expect_equal(row(2), at(25))
  everyone <- mean(c(21, 25, 33, 61, 67, 72, 78))
  study1 <- mean(c(21, 25, 33))
  imputed <- ctsem:::.ctTIImpute("age", c(21, 25, 33, NA, 61, 67, 72, 78, NA),
    NULL, includeOutcome = FALSE, group = rep(1:3, c(4, 4, 1)))
  expect_gt(imputed$mu[1L], study1)
  expect_lt(imputed$mu[1L], everyone)
  expect_equal(imputed$mu[2L], everyone)
  expect_gt(imputed$sigma[2L], imputed$sigma[1L])
  expect_equal(row(4), at(imputed$mu[1L]))
  expect_equal(row(9), at(everyone))
})

test_that("a smooth covariate is imputed from the person's other data", {
  # A covariate that is close to a time dependent predictor's mean is
  # recovered from it for the persons who lack it.
  set.seed(3)
  n <- 60
  td <- rnorm(n)
  persons <- data.frame(id = 1:n, age = 40 + 10 * td + rnorm(n, 0, 1))
  persons$age[1:5] <- NA
  long <- persons[rep(1:n, each = 3), ]
  long$time <- rep(0:2, n)
  long$Y1 <- rnorm(3 * n)
  long$TD1 <- rep(td, each = 3)
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    LAMBDA = diag(1), CINT = matrix("cint"), MANIFESTMEANS = matrix(0),
    TDpredNames = "TD1", TIpredNames = ~ s(age, k = 5)))
  out <- suppressMessages(ctsem:::.ctTIsmoothData(m, long,
    includeOutcome = FALSE))
  basis <- out$model$TIsmooths$age$basis
  imputed <- out$data$age_lin[match(1:5, out$data$id)]
  truth <- ctsem:::.ctTIsmoothEval(basis, 40 + 10 * td[1:5])[, 1]
  expect_lt(max(abs(imputed - truth)), 0.1)
})

test_that("the imputation model is each person's, from all they have", {
  set.seed(2)
  n <- 200
  x <- rnorm(n)
  z <- rnorm(n)
  target <- 2 * x - z + rnorm(n, 0, 0.3)
  truth <- target
  target[1:3] <- NA
  inputs <- cbind(x = x, z = z)
  inputs[2, "z"] <- NA
  inputs[3, ] <- NA
  out <- ctsem:::.ctTIImpute("t", target, inputs, includeOutcome = FALSE)
  expect_equal(out$person, 1:3)
  # Both inputs, one, none: each person's own model, and its own residual sd.
  expect_lt(abs(out$mu[1L] - (2 * x[1] - z[1])), 0.2)
  expect_lt(abs(out$mu[2L] - 2 * x[2]), 0.4)
  expect_equal(out$mu[3L], mean(target, na.rm = TRUE))
  expect_lt(out$sigma[1L], out$sigma[2L])
  expect_lt(out$sigma[2L], out$sigma[3L])
  expect_equal(out$sigma[3L], stats::sd(target, na.rm = TRUE))
  # A column the data cannot separate from another is dropped, not fatal.
  twin <- cbind(inputs, x2 = 2 * inputs[, "x"])
  same <- ctsem:::.ctTIImpute("t", target, twin, includeOutcome = FALSE)
  expect_equal(same$mu, out$mu)
  # Held within the observed range.
  far <- inputs
  far[1, "x"] <- 50
  expect_equal(ctsem:::.ctTIImpute("t", target, far,
    includeOutcome = FALSE)$mu[1L], max(target, na.rm = TRUE))
})

# One fit: a random intercept whose mean is a curve in age, and a model whose
# initial level varies with it, as a stationary process's does.
.tismooth_data <- function() {
  set.seed(11)
  n <- 120
  tp <- 10
  age <- runif(n, 20, 80)
  curve <- function(a) 0.6 * sin(2 * pi * (a - 20) / 60)
  drift <- -0.5
  A <- exp(drift)
  Q <- (1 - A^2) / (-2 * drift)
  rows <- lapply(seq_len(n), function(i) {
    cint <- curve(age[i]) + stats::rnorm(1, 0, 0.2)
    mu <- -cint / drift
    eta <- mu + stats::rnorm(1, 0, sqrt(1 / (-2 * drift)))
    y <- numeric(tp)
    for (t in seq_len(tp)) {
      if (t > 1) eta <- A * eta + (1 - A) * mu + stats::rnorm(1, 0, sqrt(Q))
      y[t] <- eta + stats::rnorm(1, 0, 0.4)
    }
    data.frame(id = i, time = seq_len(tp) - 1, Y1 = y, age = age[i])
  })
  list(data = do.call(rbind, rows), curve = curve)
}

.tismooth_fit <- function() fit_cached("tismooth_fit", {
  dat <- .tismooth_data()$data
  m <- ctModel(type = "ct", n.latent = 1, n.manifest = 1, LAMBDA = diag(1),
    DRIFT = matrix("drift"), DIFFUSION = matrix("diff"),
    MANIFESTVAR = matrix("mvar"), CINT = matrix("cint||||age"),
    T0MEANS = matrix("t0m||||"), T0VAR = matrix("t0v"),
    MANIFESTMEANS = matrix(0), TIpredNames = ~ s(age, k = 8), silent = TRUE)
  suppressWarnings(suppressMessages(ctFit(dat, m, backend = "julia",
    cores = 1)))
})

test_that("a curve in age is recovered, with its smoothness estimated", {
  skip_without_julia()
  fit <- .tismooth_fit()
  truth <- .tismooth_data()$curve
  expect_true(fit$optim$converged)
  expect_equal(nrow(fit$smooths), 1L)
  expect_equal(fit$smooths$param, "cint")
  # Neither a line nor every column: the curve is a sine.
  expect_gt(fit$smooths$edf, 2)
  expect_lt(fit$smooths$edf, 7)
  # The smoothing parameter is at its fixed point at the reported estimate.
  penalty <- fit$model_spec$penalty
  again <- ctsem:::.ctTIsmoothLambda(penalty, fit$uncertainty$hessian,
    as.numeric(fit$estimate$raw), as.numeric(fit$optim$gradient))
  expect_equal(again$lambda, penalty$lambda, tolerance = 0.02)
  # The fitted curve against the true one, both centred over the grid: CINT is
  # the parameter's own value at the average person.
  grid <- seq(25, 75, by = 5)
  X <- ctsem:::.ctTIsmoothEval(fit$modelbase$TIsmooths$age$basis, grid)
  fitted <- vapply(seq_along(grid), function(g) suppressMessages(as.numeric(
    ctBackendParMatrices(fit, tipreds = X[g, ])$CINT[1, 1])), 0)
  error <- (fitted - mean(fitted)) - (truth(grid) - mean(truth(grid)))
  expect_lt(sqrt(mean(error^2)), 0.15)
  expect_gt(stats::cor(fitted, truth(grid)), 0.9)
  # The basis and the smoothing parameter go back on the model, so a refit
  # evaluates the same function and starts from that smoothness.
  expect_equal(fit$modelbase$TIsmooths$age$lambda$cint, penalty$lambda)
})

# Where the finish walks the Krylov model, the exact Hessian being dear, the
# first stage ends on its walk as the smoothing passes do, and the one exact
# Hessian is the certified stage's once the smoothing parameter settles:
# certified after the first stage too, the affect pilot formed an 8-hour
# Hessian for a smoothing parameter its first pass then moved. A random drift
# keeps a Laplace coordinate, which the finish's Krylov model needs.
#
# One exact Hessian in the whole fit. Until 2026-11-15 there were three, and
# the count this test read showed two: the run dropped every walk its finish
# ended on, so the certification formed its own Hessian after the first stage
# (which `newton_hessians` does not count) and a certified stage followed each
# pass; and the pass after the certified stage's check walked and certified
# again, where it now starts on that stage's Hessian and keeps it.
test_that("a smooth fit on the Krylov finish forms no exact Hessian before its smoothing settles", {
  skip_without_julia()
  dat <- .tismooth_data()$data
  m <- ctModel(type = "ct", n.latent = 1, n.manifest = 1, LAMBDA = diag(1),
    DRIFT = matrix("drift"), DIFFUSION = matrix("diff"),
    MANIFESTVAR = matrix("mvar"), CINT = matrix("cint||||age"),
    T0MEANS = matrix("t0m||||"), T0VAR = matrix("t0v"),
    MANIFESTMEANS = matrix(0), TIpredNames = ~ s(age, k = 8), silent = TRUE)
  m$pars$indvarying <- m$pars$param %in% "drift"
  # The certification's own Hessians, formed where a run hands back none; and
  # the smoothing steps taken on the Krylov model, each with the smooth's own
  # columns exact in it. And the stages on the smoothing objective (those that
  # keep its curvature pairs): the certified one goes on from stages that
  # stopped by their own rule, so it is run as a resumed stage, whose L-BFGS
  # gives way to the finish once it stops gaining. Run as a first stage, the
  # affect pilot's crept for hours towards a hand-over 360 iterations away.
  own <- 0L
  endgame <- ctsem:::.ctBackendEndgameAt
  exact <- ctsem:::.ctTIsmoothExactColumns
  optimise <- ctsem:::.ctJuliaOptimise
  moved <- numeric()
  stages <- list()
  testthat::local_mocked_bindings(.ctBackendEndgameAt = function(...) {
    own <<- own + 1L
    endgame(...)
  }, .ctTIsmoothExactColumns = function(spec, hessian, theta, cores) {
    out <- exact(spec, hessian, theta, cores)
    cols <- ctsem:::.ctTIsmoothOwnParameters(spec, length(theta))
    moved <<- c(moved, max(abs(out - hessian)[, cols]))
    out
  }, .ctJuliaOptimise = function(...) {
    asked <- list(...)
    if (isTRUE(asked[["pairs"]])) stages[[length(stages) + 1L]] <<- list(
      certified = identical(asked[["walked"]], FALSE) || is.null(asked[["walked"]]),
      resumed = isTRUE(asked[["resumed"]]))
    optimise(...)
  }, .package = "ctsem")
  set.seed(1)
  fit <- suppressWarnings(suppressMessages(ctFit(dat, m, backend = "julia",
    intoverpop = "laplace", cores = 1,
    optimcontrol = list(checks = FALSE, newton = "krylov"))))
  expect_true(fit$optim$converged)
  expect_gt(fit$optim$newton_krylov_models, 0L)
  penalty <- fit$model_spec$penalty
  expect_gte(penalty$passes, 1L)
  expect_equal(fit$optim$newton_hessians, 1L)
  expect_equal(own, 0L)
  expect_gte(length(moved), 1L)
  expect_true(all(moved > 0))
  certified <- vapply(stages, function(s) s$certified, NA)
  resumed <- vapply(stages, function(s) s$resumed, NA)
  expect_false(certified[1L])
  expect_gte(sum(certified), 1L)
  expect_identical(resumed, certified)
  again <- ctsem:::.ctTIsmoothLambda(penalty, fit$uncertainty$hessian,
    as.numeric(fit$estimate$raw), as.numeric(fit$optim$gradient))
  expect_equal(again$lambda, penalty$lambda, tolerance = 0.02)
})

# A resume from a checkpoint that holds the certified stage's Hessian goes on
# as the run that formed it would. Where that Hessian's check moves a
# smoothing parameter, the passes walk on it at the new precisions and a
# certified stage seeded with it follows. Read as a first stage's Hessian, as it was until
# 2026-11-16, the check sent the resume through a walked pass and a second
# exact Hessian: 8 hours on the affect pilot's tied fit. The checkpoint says
# which stage wrote it (`settled`).
test_that("a resumed smooth fit takes its kept Hessian's check as the run that formed it would", {
  skip_without_julia()
  path <- tempfile(fileext = ".rds")
  side <- paste0(path, ".hessian")
  on.exit(unlink(c(path, side)), add = TRUE)
  dat <- .tismooth_data()$data
  m <- ctModel(type = "ct", n.latent = 1, n.manifest = 1, LAMBDA = diag(1),
    DRIFT = matrix("drift"), DIFFUSION = matrix("diff"),
    MANIFESTVAR = matrix("mvar"), CINT = matrix("cint||||age"),
    T0MEANS = matrix("t0m||||"), T0VAR = matrix("t0v"),
    MANIFESTMEANS = matrix(0), TIpredNames = ~ s(age, k = 8), silent = TRUE)
  m$pars$indvarying <- m$pars$param %in% "drift"
  fitto <- function() {
    set.seed(1)
    suppressWarnings(suppressMessages(ctFit(dat, m, backend = "julia",
      intoverpop = "laplace", cores = 1,
      optimcontrol = list(checks = FALSE, newton = "krylov", checkpoint = path))))
  }
  first <- fitto()
  expect_true(first$optim$converged)
  saved <- readRDS(path)
  kept <- readRDS(side)
  expect_true(isTRUE(saved$settled))
  expect_identical(kept$id, saved$hessian_id)
  # The checkpoint as it would be had the stage ended at three times the
  # smoothing parameter: the same point, and the same Hessian with the
  # penalty's precision changed on the penalised coefficients, which is the
  # exact one there. Its check moves the parameter back.
  penalty <- first$model_spec$penalty
  expect_equal(saved$state$lambda, penalty$lambda)
  index <- as.integer(unlist(penalty$index))
  more <- rep(2 * penalty$lambda, lengths(penalty$index))
  saved$state$lambda <- 3 * saved$state$lambda
  kept$state <- saved$state
  kept$hessian[cbind(index, index)] <- kept$hessian[cbind(index, index)] - more
  saveRDS(saved, path)
  saveRDS(kept, side)
  stages <- list()
  optimise <- ctsem:::.ctJuliaOptimise
  testthat::local_mocked_bindings(.ctJuliaOptimise = function(...) {
    asked <- list(...)
    stages[[length(stages) + 1L]] <<- list(
      walked = !(identical(asked[["walked"]], FALSE) || is.null(asked[["walked"]])),
      held = !is.null(asked[["hessian"]]))
    optimise(...)
  }, .package = "ctsem")
  second <- fitto()
  expect_true(second$optim$converged)
  # The first stage on the kept Hessian; then the check's move as passes
  # that walk on it until the parameters settle, and one certified stage
  # handed it at the settled precisions. Each holds a Hessian, and the
  # first and last do not end on a walk.
  expect_gte(length(stages), 3L)
  walked <- vapply(stages, function(s) s$walked, NA)
  expect_false(walked[1L])
  expect_false(walked[length(walked)])
  expect_true(any(walked))
  expect_true(all(vapply(stages, function(s) s$held, NA)))
  expect_gte(second$model_spec$penalty$passes, 1L)
  expect_equal(second$model_spec$penalty$lambda, penalty$lambda, tolerance = 0.05)
  expect_equal(second$estimate$loglik, first$estimate$loglik, tolerance = 1e-3)
})

test_that("summary() and ctTIpredEffects() report a smooth whole", {
  skip_without_julia()
  fit <- .tismooth_fit()
  s <- suppressMessages(summary(fit))
  expect_equal(rownames(s$smooths), "s(age) on cint")
  expect_lt(s$smooths$p, 0.01)
  # Most of the person differences in cint are the curve's: a sine of
  # amplitude 0.6 against a random sd of 0.2.
  expect_gt(s$smooths$share, 0.5)
  expect_lt(s$smooths$share, 1)
  expect_false(any(grepl("^tip_age_", rownames(s$tipreds))))
  expect_output(print(s), "Smooth time-independent predictor terms")
  # AIC counts the smooth by its effective degrees of freedom.
  edf <- fit$model_spec$penalty$edf
  expect_equal(s$aic, 2 * (s$npars - fit$model_spec$penalty$m + edf) -
    2 * fit$estimate$loglik, tolerance = 1e-6)
  effects <- suppressMessages(ctTIpredEffects(fit, whichTIpreds = "age",
    parmatrices = FALSE, whichpars = "all", nsubjects = 15))
  expect_equal(colnames(effects$x), "age")
  expect_equal(dim(effects$y)[1L], 15L)
  expect_true("cint" %in% dimnames(effects$y)[[2L]])
  # The bands carry the coefficients' uncertainty.
  width <- effects$y[, "cint", 3L] - effects$y[, "cint", 1L]
  expect_true(all(width > 0))
  # Each person as a point at their own age.
  persons <- effects$persons[effects$persons$variable == "cint", ]
  expect_equal(nrow(persons), 120L)
  g <- suppressMessages(ctTIpredEffects(fit, whichTIpreds = "age",
    parmatrices = FALSE, nsubjects = 15, plot = TRUE))
  expect_s3_class(g, "ggplot")
})

test_that("a fit given as inits starts a changed model where it stopped", {
  skip_without_julia()
  fit <- .tismooth_fit()
  # The same model without the smooth: its parameters are the fit's, by name.
  dat <- .tismooth_data()$data
  plain <- ctModel(type = "ct", n.latent = 1, n.manifest = 1, LAMBDA = diag(1),
    DRIFT = matrix("drift"), DIFFUSION = matrix("diff"),
    MANIFESTVAR = matrix("mvar"), CINT = matrix("cint"),
    T0MEANS = matrix("t0m"), T0VAR = matrix("t0v"),
    MANIFESTMEANS = matrix(0), silent = TRUE)
  spec <- suppressMessages(ctFit(dat, plain, backend = "julia", fit = FALSE))
  npar <- ctsem:::.ctBackendNpar(spec)
  start <- expect_message(ctsem:::.ctJuliaInitsFromFit(fit, spec, npar),
    paste0("for ", npar, " of ", npar))
  names <- ctsem:::.ctBackendRawParameterNames(list(model_spec = spec), npar)
  old <- ctsem:::.ctBackendRawParameterNames(fit, length(fit$estimate$raw))
  expect_equal(start, as.numeric(fit$estimate$raw)[match(names, old)])
  # And the smooth model from the plain one's values: its coefficients at zero.
  back <- suppressMessages(ctsem:::.ctJuliaInitsFromFit(
    structure(list(model_spec = spec, estimate = list(raw = start)),
      class = c("ctJuliaFit", "ctFit")), fit$model_spec,
    length(fit$estimate$raw)))
  coefficients <- unlist(fit$model_spec$penalty$index)
  expect_true(all(back[coefficients] == 0))
  expect_error(suppressMessages(ctFit(dat, plain, backend = "stan",
    inits = fit, fit = FALSE)), "backend='julia' fit")
})

test_that("the Laplace Hessian takes a smooth's columns from its target's", {
  skip_without_julia()
  fit <- .tismooth_fit()
  dat <- .tismooth_data()$data
  spec <- suppressMessages(ctFit(dat, fit$modelbase, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  npar <- ctsem:::.ctBackendNpar(spec)
  x <- suppressMessages(ctsem:::.ctJuliaInitsFromFit(fit, spec, npar))
  objective <- ctsem:::.ctJuliaObjective(structure(spec,
    class = c("ctJuliaModel", "ctFitModel")), hybrid = TRUE)
  module <- ctsem:::.ctJuliaModule(spec$project)
  ctsem:::.ctJuliaEval(paste0("ctsem_test_nderived(o, n) = ",
    "length(ContinuousTimeSEM._laplace_ti_derivable(o, collect(1:n)))"))
  derived <- ctsem:::.ctJuliaCall("ctsem_test_nderived", objective,
    as.integer(npar))
  # Every coefficient of the smooth: seven, on cint.
  expect_equal(as.integer(derived), 7L)
  fast <- ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_hessian(objective,
    ctsem:::.ctJuliaNumericVector(x)))
  slow <- ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_hessian(objective,
    ctsem:::.ctJuliaNumericVector(x), localised = FALSE))
  expect_equal(as.matrix(fast), as.matrix(slow), tolerance = 1e-5)
})

test_that("smoothing parameters are settled when the curve's degrees of freedom are", {
  # A term at its linear limit: its coefficients at zero. From anywhere
  # beyond the limit its lambda runs to its bound at no change in the curve,
  # which the effective degrees of freedom say and a change in lambda does
  # not. No fit.
  penalty <- data.frame(lambda = 1e3)
  penalty$index <- list(1:3)
  info <- diag(c(1, 1, 1, 5e5))
  theta <- c(0, 0, 0, 0.3)
  gfs <- ctsem:::.ctTIsmoothLambda(penalty,
    -(info + diag(c(rep(penalty$lambda, 3), 0))), theta)
  expect_gt(abs(log(gfs$lambda / penalty$lambda)), 1)
  expect_lt(ctsem:::.ctTIsmoothMoved(gfs), 0.01)
  # The bound is the sharpest curvature of the rest of the model, so the
  # penalised information's largest eigenvalue stays the data's to a factor
  # of two, and a rule relative to it reads what it did without the smooth.
  expect_equal(gfs$lambda, 5e5)
  expect_lte(max(eigen(-gfs$hessian, symmetric = TRUE, only.values = TRUE)$values),
    2 * max(diag(info)))
  # The same term from a weak penalty has its three degrees of freedom to lose.
  penalty$lambda <- 1
  gfs <- ctsem:::.ctTIsmoothLambda(penalty,
    -(info + diag(c(rep(penalty$lambda, 3), 0))), theta)
  expect_gt(ctsem:::.ctTIsmoothMoved(gfs), 1)
})

test_that("a curvature model takes the smooths' own columns exactly", {
  # No fit: the columns' source is mocked. One term on parameter 1: its
  # penalised coefficients (2, 3), its line (4) and another predictor's
  # coefficient on the same parameter (5) are its own, with the parameter;
  # 6 and 7 are the rest, 7 with a coefficient of its own that is not.
  penalty <- data.frame(lambda = 5, linear = 4L, parameter = 1L)
  penalty$index <- list(2:3)
  spec <- list(penalty = penalty, ti_effects = data.frame(
    coefficient = c(2L, 3L, 4L, 5L, 7L), parameter = c(1L, 1L, 1L, 1L, 6L)))
  own <- 1:5
  rest <- 6:7
  expect_identical(ctsem:::.ctTIsmoothOwnParameters(spec, 7L), own)
  truth <- -(crossprod(matrix(sin(1:49), 7, 7)) + diag(7))
  model <- -(crossprod(matrix(cos(1:49), 7, 7)) + diag(7))
  asked <- NULL
  testthat::local_mocked_bindings(.ctBackendHessianColumns = function(spec, est, columns) {
    asked <<- columns
    truth[, columns, drop = FALSE]
  }, .ctBackendWithMaxChunks = function(chunks, expr) expr, .package = "ctsem")
  out <- ctsem:::.ctTIsmoothExactColumns(spec, model, numeric(7))
  expect_identical(asked, own)
  # The own columns are the exact ones, the matrix symmetric and a maximum's,
  # and the rest keeps the covariance the model gave it.
  expect_equal(out[, own], truth[, own])
  expect_equal(out, t(out))
  expect_true(all(eigen(-out, symmetric = TRUE, only.values = TRUE)$values > 0))
  expect_equal(solve(-out)[rest, rest], solve(-model)[rest, rest])
  expect_gt(max(abs(out[rest, rest] - model[rest, rest])), 0.01)
  # Columns that did not solve, or an own block that is no maximum's, leave
  # the model as it came.
  testthat::local_mocked_bindings(.ctBackendHessianColumns = function(spec, est, columns)
    matrix(NaN, 7, length(columns)), .package = "ctsem")
  expect_identical(ctsem:::.ctTIsmoothExactColumns(spec, model, numeric(7)), model)
  testthat::local_mocked_bindings(.ctBackendHessianColumns = function(spec, est, columns)
    -truth[, columns, drop = FALSE], .package = "ctsem")
  expect_identical(ctsem:::.ctTIsmoothExactColumns(spec, model, numeric(7)), model)
})
