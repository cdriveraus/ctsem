# ctLOO's two approximations to leaving data out on julia fits: PSIS
# (method = 'psis') on the filter's route and on sampled fits, and
# casewiseApproximation on every fold type. Each against an answer computed
# another way: a grid over the two free population parameters for PSIS's
# Bayesian leave-one-out, refitting for casewise, and finite differences for
# the score casewise steps by. The Laplace route's PSIS is checked the same way
# in test-julia-loo-laplace.R.

skip_on_cran()

# A linear Gaussian random intercept with every dynamic parameter fixed, so
# only the intercept's mean and population sd are free: the exact
# leave-one-subject-out density is a two-dimensional integral. Subject 1's
# intercept is far out, so it is the unit that moves the estimate.
.lapx_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), DRIFT = matrix(-0.5), DIFFUSION = matrix(0.6),
    MANIFESTVAR = matrix(0.4), MANIFESTMEANS = matrix("mmean"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE
  model
}

# Simulated here rather than through ctGenerate, whose draw stream moves under
# unrelated commits.
.lapx_data <- function(nsubjects = 16, nobs = 5, seed = 11) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(nsubjects), function(i) {
    intercept <- if (i == 1) 6 else stats::rnorm(1, 1, 1.2)
    state <- stats::rnorm(1, 0, 0.5)
    y <- numeric(nobs)
    for (t in seq_len(nobs)) {
      if (t > 1) {
        state <- exp(-0.5) * state + stats::rnorm(1, 0, sqrt(0.36 * (1 - exp(-1))))
      }
      y[t] <- state + intercept + stats::rnorm(1, 0, 0.4)
    }
    data.frame(id = i, time = seq_len(nobs) - 1, Y1 = y)
  }))
}

.lapx_cache <- new.env(parent = emptyenv())
.lapx_get <- function(name, make) {
  if (is.null(.lapx_cache[[name]])) .lapx_cache[[name]] <- make()
  .lapx_cache[[name]]
}
.lapx_fit <- function() .lapx_get("fit", function()
  suppressMessages(ctFit(.lapx_data(), .lapx_model(), backend = "julia",
    intoverpop = "augmented", priors = TRUE, verbose = 0)))
.lapx_sampled <- function() .lapx_get("sampled", function()
  suppressWarnings(suppressMessages(ctFit(.lapx_data(), .lapx_model(),
    backend = "julia", intoverpop = "augmented", priors = TRUE,
    optimize = FALSE, verbose = 0))))

# Exact: log p(y_i | y_-i) = log int p(y | th) dth - log int p(y_-i | th) dth,
# on a grid in the covariance's own axes wide enough that the tails beyond it
# do not matter. The filter's subject terms are exact marginals for this
# model, and the log prior is the objective's own, so this is the target PSIS
# estimates, from an optimized fit or a sampled one.
.lapx_exact <- function() .lapx_get("exact", function() {
  fit <- .lapx_fit()
  eig <- eigen(fit$estimate$cov, symmetric = TRUE)
  axis <- seq(-9, 9, length.out = 91)
  nodes <- as.matrix(expand.grid(a = axis, b = axis))
  grid <- as.numeric(fit$estimate$raw) + eig$vectors %*% (sqrt(eig$values) * t(nodes))
  terms <- ctsem:::.ctBackendUnitTerms(ctsem:::.ctBackendSpec(fit), grid)
  stopifnot(all(is.finite(terms$value)))
  lse <- ctsem:::.ctLogSumExp
  vapply(seq_len(nrow(terms$unit_loglik)), function(i)
    lse(terms$value) - lse(terms$value - terms$unit_loglik[i, ]), numeric(1))
})

rms <- function(x) sqrt(mean(x^2))

test_that("PSIS on the filter's route matches the exact integral over the population", {
  skip_without_julia()
  skip_if_not_installed("loo")
  fit <- .lapx_fit()
  expect_false(.ctBackendIsLaplace(fit))
  expect_length(fit$estimate$raw, 2L)
  exact <- .lapx_exact()

  set.seed(1)
  psis <- suppressWarnings(ctLOO(fit, method = "psis", nsamples = 1000))
  set.seed(1)
  more <- suppressWarnings(ctLOO(fit, method = "psis", nsamples = 4000))
  expect_s3_class(psis, "ctLOOpsis")
  expect_identical(psis$level, "id")
  expect_equal(nrow(psis$pointwise), 16L)
  expect_equal(psis$ndraws, 1000L)
  # Absolute, at two draw counts: a wrong estimator can sit close to the answer
  # at one count, but its error does not shrink with more draws. Locally the
  # rms was 0.010 at 1000 draws and 0.007 at 4000, the worst unit 0.024.
  expect_lt(rms(psis$pointwise$elpd_loo - exact), 0.025)
  expect_lt(rms(more$pointwise$elpd_loo - exact), 0.015)
  expect_lt(max(abs(more$pointwise$elpd_loo - exact)), 0.05)
  expect_true(all(more$pointwise$pareto_k < 0.7))
  expect_false(any(more$pointwise$refined))
  # The outlier is the unit whose prediction suffers most from leaving it out.
  expect_identical(which.max(more$pointwise$p_loo), 1L)
})

test_that("a unit redrawn around its leave-one-out step matches the exact integral", {
  skip_without_julia()
  skip_if_not_installed("loo")
  fit <- .lapx_fit()
  exact <- .lapx_exact()
  set.seed(2)
  again <- ctsem:::.ctPsisRefine(fit, ctsem:::.ctBackendSpec(fit), 1:3,
    as.numeric(fit$estimate$raw), fit$estimate$cov, 2000L, 1.5)
  expect_true(all(is.finite(again$elpd)))
  expect_true(all(again$k < 0.7))
  expect_lt(max(abs(again$elpd - exact[1:3])), 0.05)
})

test_that("PSIS on a sampled fit uses its posterior draws", {
  skip_without_julia()
  skip_if_not_installed("loo")
  sampled <- .lapx_sampled()
  expect_identical(sampled$sample$target, "marginal")
  exact <- .lapx_exact()
  psis <- suppressWarnings(ctLOO(sampled, method = "psis"))
  expect_equal(psis$ndraws, nrow(sampled$estimate$rawposterior))
  expect_true(is.na(psis$proposal_k))
  expect_match(psis$note, "from the posterior draws")
  # Monte Carlo error of a few hundred MCMC draws; locally 0.021.
  expect_lt(rms(psis$pointwise$elpd_loo - exact), 0.05)
  thinned <- suppressWarnings(ctLOO(sampled, method = "psis", nsamples = 200))
  expect_lte(thinned$ndraws, 200L)
  # And the row level, which conditions on one parameter vector, is refused.
  expect_error(ctLOO(sampled, method = "psis", subjectwise = FALSE),
    "leaves one row out only on an optimized fit")
})

test_that("casewise steps follow refits on row and subject folds, and a bar refits the far ones", {
  skip_without_julia()
  fit <- .lapx_fit()
  rowsubject <- .ctFitRowSubject(fit)

  # Row folds: each fold leaves its subjects partly observed, so the step is
  # the difference of the full and training data's scores.
  set.seed(3)
  stepped <- suppressMessages(ctLOO(fit, folds = 5, subjectwise = FALSE,
    casewiseApproximation = TRUE, cores = 1))
  set.seed(3)
  refit <- suppressMessages(ctLOO(fit, folds = 5, subjectwise = FALSE, cores = 1))
  expect_identical(stepped$foldrows, refit$foldrows)
  expect_length(stepped$foldinfluence, 5L)
  expect_true(all(is.finite(stepped$foldinfluence) & stepped$foldinfluence > 0))
  expect_true(all(stepped$foldstepped))
  expect_match(stepped$scoring, "one Newton step")
  # Locally the largest row difference was 0.014.
  expect_lt(max(abs(stepped$outsampleLogLikRow - refit$outsampleLogLikRow),
    na.rm = TRUE), 0.05)
  # Stepped folds are scored together, each subject at the parameters of the
  # fold holding its rows -- several passes here, as a subject's rows span
  # folds: the held-out rows are what a full pass at the fold's parameters
  # gives them, and no other row is scored.
  heldrows <- function(loo, k) {
    rows <- loo$foldrows[[k]]
    full <- .ctBackendRowLoglik(fit, loo$foldpars[, k])
    full[full == 0] <- NA
    expect_equal(loo$LogLikRowFolds[[k]][rows], full[rows], tolerance = 1e-10)
    expect_true(all(is.na(loo$LogLikRowFolds[[k]][-rows])))
  }
  for (k in seq_along(stepped$foldrows)) heldrows(stepped, k)

  # A bar steps the folds whose influence is under it and refits the rest,
  # which then agree with the plain refits to the optimiser's tolerance.
  bar <- stats::median(stepped$foldinfluence)
  set.seed(3)
  mixed <- suppressMessages(ctLOO(fit, folds = 5, subjectwise = FALSE,
    casewiseApproximation = bar, cores = 1))
  expect_identical(mixed$foldstepped, stepped$foldinfluence < bar)
  for (k in which(!mixed$foldstepped)) {
    rows <- mixed$foldrows[[k]]
    expect_lt(max(abs(mixed$outsampleLogLikRow[rows] - refit$outsampleLogLikRow[rows]),
      na.rm = TRUE), 5e-3)
  }

  # Subject folds: the outlier moves the estimate furthest, and that is where
  # one step falls short of the refit; elsewhere it is close.
  set.seed(4)
  stepped <- suppressMessages(ctLOO(fit, folds = 16, casewiseApproximation = TRUE,
    cores = 1))
  set.seed(4)
  refit <- suppressMessages(ctLOO(fit, folds = 16, cores = 1))
  # One pass for all sixteen folds here.
  for (k in c(1L, 9L, 16L)) heldrows(stepped, k)
  bysubject <- function(x) as.numeric(tapply(x$outsampleLogLikRow, rowsubject, sum))
  error <- abs(bysubject(stepped) - bysubject(refit))
  foldsubject <- vapply(stepped$foldrows, function(r) rowsubject[r[1]], numeric(1))
  influence <- stepped$foldinfluence[order(foldsubject)]
  # The error grows with the influence (locally 0.46 at 0.97, 0.04 at 0.23,
  # at most 0.01 below 0.11), which is what makes the influence a bar.
  expect_identical(which.max(influence), 1L)
  expect_identical(which.max(error), 1L)
  expect_gt(stats::cor(influence, error), 0.8)
  expect_lt(max(error[influence < 0.1]), 0.02)
})

# Subjects in studies: a subject fold leaves its study partly observed, so the
# step needs the training data's score for that study.
.lapx_nested_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), MANIFESTMEANS = matrix("mmean"),
    id = c("subject", "study"))))
  model$pars$indvarying <- model$pars$param %in% "mmean"
  model
}
.lapx_nested_data <- function(nstudy = 6, npersub = 4, nobs = 5) {
  set.seed(20260825)
  rows <- list(); sid <- 0
  for (g in seq_len(nstudy)) {
    studyeffect <- stats::rnorm(1, 0, 0.3)
    for (j in seq_len(npersub)) {
      sid <- sid + 1
      intercept <- 1 + studyeffect + stats::rnorm(1, 0, 0.5)
      state <- stats::rnorm(1, 0, 0.5); out <- numeric(nobs)
      for (t in seq_len(nobs)) {
        if (t > 1) {
          state <- exp(-0.4) * state +
            stats::rnorm(1, 0, sqrt(0.36 / 0.8 * (1 - exp(-0.8))))
        }
        out[t] <- state + intercept + stats::rnorm(1, 0, 0.3)
      }
      rows[[length(rows) + 1L]] <- data.frame(subject = sid, study = g,
        time = seq_len(nobs) - 1, Y1 = out)
    }
  }
  do.call(rbind, rows)
}

test_that("casewise inside a group steps by the subject's share of its group's score", {
  skip_without_julia()
  fit <- suppressWarnings(suppressMessages(ctFit(.lapx_nested_data(),
    .lapx_nested_model(), backend = "julia", intoverpop = "laplace",
    priors = TRUE, verbose = 0)))
  spec <- .ctBackendSpec(fit)
  est <- .ctLaplaceOptimum(fit)
  units <- .ctBackendUnitTerms(spec, est)
  expect_equal(nrow(units$unit_loglik), 6L)

  # What leaving subject 3 out removes is its contribution to its study's
  # marginal, the full study's term less the study's without it; the step's
  # score is that contribution's gradient.
  rowsubject <- .ctFitRowSubject(fit)
  heldout <- which(rowsubject == 3)
  training <- spec$data
  training[heldout, "Y1"] <- NA
  trainfit <- .ctFitReplaceData(fit, training)
  ctx <- list(scores = .ctBackendScoreMatrix(fit, est, prior = FALSE),
    covariance = fit$estimate$cov, unitOf = units$unit, rowsubject = rowsubject,
    est = est)
  removed <- .ctLOOFoldScore(ctx, heldout, training = function() trainfit)
  U <- units$unit[3]
  contribution <- function(p) {
    .ctBackendUnitTerms(spec, p)$unit_loglik[U, 1] -
      .ctBackendUnitTerms(trainfit$model_spec, p)$unit_loglik[U, 1]
  }
  h <- 1e-4
  fd <- vapply(seq_along(est), function(j) {
    e <- replace(numeric(length(est)), j, h)
    (contribution(est + e) - contribution(est - e)) / (2 * h)
  }, numeric(1))
  expect_equal(removed, fd, tolerance = 1e-5)
  # Without the training function the partly observed study cannot be
  # stepped, and the fold is refit instead.
  expect_null(.ctLOOFoldScore(ctx, heldout))

  set.seed(7)
  stepped <- suppressWarnings(suppressMessages(ctLOO(fit, folds = 6,
    casewiseApproximation = TRUE, cores = 1)))
  set.seed(7)
  refit <- suppressWarnings(suppressMessages(ctLOO(fit, folds = 6, cores = 1)))
  expect_true(all(stepped$foldstepped))
  # Each fold holds out a sixth of the data, so the steps are long (locally
  # influence 0.33-1.08) and capture three quarters of the leave-one-out
  # penalty; the sign and the order of magnitude are what is pinned.
  penalty <- refit$insampleLogLik - refit$outsampleLogLik
  captured <- refit$insampleLogLik - stepped$outsampleLogLik
  expect_gt(captured, 0.5 * penalty)
  expect_lt(captured, penalty + 0.5)
})
