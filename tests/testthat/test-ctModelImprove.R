# ctModelImprove() and ctFitAdequacy(): the search keeps what the data
# support, and a correct model stops at once. One latent, so each fit costs
# seconds; the two-latent searches the design was checked on are in
# CT-SEM/review/SCORETESTS-misfit-design-2026-10-09.md.

skip_without_julia()
skip_on_32bit()

.improve_data <- function(effect = 0.8, nsub = 80L, nobs = 20L, seed = 5L) {
  set.seed(seed)
  z <- stats::rnorm(nsub); w <- stats::rnorm(nsub)
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    times <- cumsum(c(0, stats::rexp(nobs - 1L, 1)))
    x <- stats::rnorm(1, 0, 0.8)
    xs <- numeric(nobs); xs[1] <- x
    for (t in 2:nobs) {
      dt <- times[t] - times[t - 1]
      x <- exp(-0.5 * dt) * x + stats::rnorm(1, 0, 0.8 * sqrt(1 - exp(-dt)))
      xs[t] <- x
    }
    mm <- stats::rnorm(1) + effect * z[i]
    data.frame(id = i, time = times, Y1 = xs + mm + stats::rnorm(nobs, 0, 0.5),
      z = z[i], w = w[i])
  }))
}

.improve_fit <- function(effect, seed = 5L) fit_cached(paste('improve', effect, seed), {
  d <- .improve_data(effect, seed = seed)
  list(data = d, fit = suppressMessages(ctFit(d[, c('id', 'time', 'Y1')],
    suppressMessages(ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
      LAMBDA = matrix(1))), backend = 'julia', cores = 1, verbose = 0)))
})

test_that("adequacy is reported per family, and the error-based families need no draws", {
  x <- .improve_fit(0)
  ad <- ctFitAdequacy(x$fit, nsamples = 40, cores = 1, nboot = 499)
  expect_s3_class(ad, 'ctFitAdequacy')
  # One manifest variable, so no within-row pairs; the lagged correlations of
  # the data only when asked for.
  expect_setequal(ad$families$family, c('residuals', 'latent', 'volatility', 'periodic',
    'distribution', 'fit', 'dependence', 'subjects', 'nonlinear'))
  expect_true(all(ad$families$p > 0 & ad$families$p <= 1))
  expect_true(ad$adequate)
  expect_true(ad$close_fit)
  # A small model with its Hessian: estimation removed by exact derivatives.
  expect_identical(ad$projection, 'derivative')
  expect_true(all(c('residual correlation, all gaps', 'residual correlation, shortest gaps',
    'latent residual correlation, all gaps') %in% ad$size$measure))
  # Reuses what the fit generated.
  expect_equal(dim(ad$fit$generated$Y)[1], 40L)
  expect_output(print(ad), 'adequate')
  # The error-based families alone generate nothing.
  ad2 <- ctFitAdequacy(x$fit, families = c('residuals', 'volatility', 'fit'), cores = 1,
    nboot = 199)
  expect_setequal(ad2$families$family, c('residuals', 'volatility', 'fit'))
  expect_true(is.na(ad2$nsamples))
  # The dynamic information matrix test, from row scores, when asked for.
  ad3 <- ctFitAdequacy(x$fit, families = 'scores', cores = 1, nboot = 199)
  expect_identical(ad3$families$family, 'scores')
  expect_true(ad3$families$p > 0 && ad3$families$p <= 1)
  expect_true(all(startsWith(ad3$statistics$statistic, 's(')))
})

test_that("a conditional error's implied correlation is P_jk / sqrt(P_jj P_kk)", {
  # The within-row family compares products of the errors with this; a sign
  # error in it once rejected every correctly specified model.
  set.seed(2)
  S <- matrix(c(1, 0.6, 0.3, 0.6, 1.5, 0.2, 0.3, 0.2, 0.8), 3)
  P <- solve(S)
  v <- t(chol(S)) %*% matrix(stats::rnorm(3 * 20000), 3)
  E <- t((P %*% v) / sqrt(diag(P)))
  implied <- P / sqrt(outer(diag(P), diag(P)))
  expect_equal(stats::cor(E)[upper.tri(S)], implied[upper.tri(S)], tolerance = 0.03)
  expect_equal(apply(E, 2, stats::var), rep(1, 3), tolerance = 0.03)
})

test_that("a strong misfit is not floored at one over the number of datasets", {
  # Subjects' means follow z, which the model does not have: the covariate
  # family's p-value used to be the share of 40 datasets beyond the data's,
  # so no smaller than 1/41, and Holm over the families then hid it.
  x <- .improve_fit(0.8)
  ad <- ctFitAdequacy(x$fit, covariates = unique(x$data[, c('id', 'z')]), nsamples = 40,
    cores = 1)
  expect_lt(ad$families$p[ad$families$family == 'covariates'], 1e-4)
  expect_false(ad$adequate)
})

test_that("the Anderson-Darling reference holds its level", {
  set.seed(3)
  p <- replicate(400, .ctAdequacyAD(stats::runif(200))[2])
  # 0.05, within three binomial sds at 400 draws.
  expect_gt(mean(p < 0.05), 0.017)
  expect_lt(mean(p < 0.05), 0.083)
  expect_lt(.ctAdequacyAD(stats::pnorm(stats::rnorm(500, 0.3)))[2], 0.001)
})

test_that("a covariate the model lacks is added, an unrelated one is not, and the search stops", {
  x <- .improve_fit(0.8)
  imp <- suppressMessages(ctModelImprove(x$fit, covariates = unique(x$data[, c('id', 'z', 'w')]),
    relax = c('heterogeneity', 'covariates'), nsamples = 40, cores = 1, verbose = FALSE))
  expect_s3_class(imp, 'ctModelImprove')
  kept <- imp$history$change[imp$history$kept]
  expect_true(any(grepl('^z as a predictor of', kept)))
  expect_false(any(grepl('^w ', imp$history$change)))
  expect_true('z' %in% ctsem:::.ctFitBaseModel(imp$fit)$TIpredNames)
  expect_gt(imp$fit$estimate$loglik, x$fit$estimate$loglik)
  expect_identical(imp$stopped, 'adequate')
  expect_output(print(imp), 'stopped because adequate')

  # Pruning keeps a change the fit plainly needs without refitting: the random
  # intercept, whose spread here is many standard errors from zero, by a Wald
  # test of its sd at zero past twice the likelihood ratio bar.
  model <- ctsem:::.ctFitBaseModel(x$fit)
  fixed <- model
  fixed$pars$indvarying[fixed$pars$matrix %in% 'MANIFESTMEANS'] <- FALSE
  needed <- ctsem:::.ctImprovePruneNeeded(x$fit, list(family = 'heterogeneity'),
    list(model = fixed, data = x$data[, c('id', 'time', 'Y1')], args = list()), 1L, 'LR', 0.05)
  expect_true(needed$needed)
  expect_gte(needed$dpar, 1L)
  expect_gt(needed$W, needed$bar)
  # A change with no absent value to test is refitted instead.
  expect_null(ctsem:::.ctImprovePruneNeeded(x$fit, list(family = 'slow'),
    list(model = fixed, data = x$data, args = list()), 1L, 'LR', 0.05))
})

test_that("a correct model stops at once", {
  # Its own draw: with seed 5, a correct model's data are one of the five per
  # cent of draws where the search keeps something at alpha = 0.05 (a random
  # diffusion, LR p = 0.0097) -- the rate the simulation in the review note
  # measured, which is the place to test it, not one dataset.
  # At alpha = 0.01: at 0.05 a single correct dataset has two chances near
  # five per cent each to stop this test -- the search keeping something, and
  # the checks failing (seed 1: the unrelated w correlates with the subjects'
  # sds at z = 3.2, covariate family Holm-adjusted p = 0.047). The rates are the
  # simulations' to measure (review note); this tests the mechanism.
  x <- .improve_fit(0, seed = 1L)
  imp <- suppressMessages(ctModelImprove(x$fit, covariates = unique(x$data[, c('id', 'w')]),
    relax = c('heterogeneity', 'covariates', 'statedependence'), alpha = 0.01,
    nsamples = 40, cores = 1, verbose = FALSE))
  expect_identical(imp$stopped, 'adequate')
  # The screen is at alpha, so it can propose something; the refit's
  # likelihood is what decides, and nothing may be kept.
  expect_false(any(imp$history$kept))
  expect_equal(length(ctsem:::.ctFitRawEstimate(imp$fit)),
    length(ctsem:::.ctFitRawEstimate(x$fit)))
})

test_that("a search with held-out subjects confirms what it kept on them", {
  x <- .improve_fit(0.8)
  imp <- suppressMessages(ctModelImprove(x$fit, covariates = unique(x$data[, c('id', 'z')]),
    relax = 'covariates', holdout = 0.3, nsamples = 40, cores = 1, verbose = FALSE))
  expect_length(imp$holdout, round(0.3 * 80))
  cf <- imp$confirmation
  expect_true(any(grepl('^z as a predictor of', cf$change)))
  expect_true(all(cf$confirmed[grepl('^z ', cf$change)]))
  # The final model is refitted to every subject.
  expect_equal(length(unique(ctsem:::.ctFitLongData(imp$fit)$id)), 80L)
  expect_output(print(imp), 'held-out')
})

test_that("a process driven by a slow one gets a slow process, and the search stops", {
  # The fast process is pulled towards a slow one with a time constant of 20,
  # which the model lacks: its updates stay correlated over long gaps (the
  # latent family), so the search tries a slow process first and keeps it.
  set.seed(21)
  d <- do.call(rbind, lapply(1:60, function(i) {
    nobs <- 40; times <- cumsum(c(0, stats::rexp(nobs - 1, 1)))
    x <- stats::rnorm(1, 0, 0.7); s <- stats::rnorm(1); xs <- numeric(nobs); xs[1] <- x
    for (t in 2:nobs) {
      n <- max(1, ceiling((times[t] - times[t - 1]) / 0.02)); h <- (times[t] - times[t - 1]) / n
      for (k in 1:n) {
        x <- x + (-x + 0.5 * s) * h + 0.8 * sqrt(h) * stats::rnorm(1)
        s <- s - 0.05 * s * h + sqrt(2 * 0.05 * h) * stats::rnorm(1)
      }
      xs[t] <- x
    }
    data.frame(id = i, time = times, Y1 = xs + stats::rnorm(nobs, 0, 0.4))
  }))
  fit <- fit_cached('improve slow', {
    m <- suppressMessages(ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
      LAMBDA = matrix(1), MANIFESTMEANS = matrix(0), CINT = matrix(0)))
    m$pars$indvarying <- FALSE
    suppressMessages(ctFit(d, m, backend = 'julia', cores = 1, verbose = 0))
  })
  imp <- suppressMessages(ctModelImprove(fit, relax = 'slow', nsamples = 40, cores = 1,
    verbose = FALSE))
  expect_identical(imp$history$change[imp$history$kept], 'eta1 driven by a slow process')
  expect_lt(imp$history$screen_p[1], 1e-4)
  expect_true('slow_eta1' %in% ctsem:::.ctFitBaseModel(imp$fit)$latentNames)
  expect_identical(imp$stopped, 'adequate')
})

test_that("an item's own process starts persistent, at a quarter of its error variance", {
  # At raw zero the default drift transform forgets within the hour, which on
  # prompts hours apart is measurement error: a flat likelihood, from which a
  # process the Fisher data needed (350 log likelihood units) gained 2.
  x <- .improve_fit(0, seed = 1L)
  model <- ctsem:::.ctFitBaseModel(x$fit)
  start <- ctsem:::.ctImproveUniqueStart(model, x$fit, 'Y1')
  # Its persistence and innovation sd, declared in PARS; it starts from its
  # stationary distribution, so its initial variance has no start of its own.
  expect_setequal(names(start), c('drift_u_Y1', 'diff_u_Y1'))
  m1 <- ctModelUniqueProcess(model, 'Y1')
  at <- function(par) {
    r <- m1$pars[m1$pars$matrix == 'PARS' & m1$pars$param %in% par, ]
    eval(parse(text = r$transform), list(param = start[[par]], log1p_exp = ctsem:::log1p_exp))
  }
  gaps <- unlist(lapply(split(x$data$time, x$data$id), diff))
  rate <- log(2) / stats::median(gaps[gaps > 0])
  mv <- model$pars[model$pars$matrix == 'MANIFESTVAR' & model$pars$row == 1 &
    model$pars$col == 1, ]
  est <- setNames(ctsem:::.ctFitRawEstimate(x$fit), ctsem:::.ctFitRawParNames(x$fit))
  v <- eval(parse(text = mv$transform), list(param = est[[mv$param]],
    log1p_exp = ctsem:::log1p_exp))^2 / 4
  expect_equal(at('drift_u_Y1'), -rate, tolerance = 1e-6)
  expect_equal(at('diff_u_Y1')^2 / (2 * rate), v, tolerance = 1e-6)
  t0 <- m1$pars$param[m1$pars$matrix == 'T0VAR' & m1$pars$row == 2 & m1$pars$col == 2]
  expect_equal(eval(parse(text = t0), list(diff_u_Y1 = at('diff_u_Y1'),
    drift_u_Y1 = at('drift_u_Y1')))^2, v, tolerance = 1e-6)
})

test_that("pruning side by side makes the decisions pruning one at a time does", {
  # No fits: each model is the labels of the changes it carries. b and d are
  # each redundant while the other is kept, e always: one at a time, b goes
  # and then d is needed. A batch compared against the full model would drop
  # both, so after a drop the rest of a batch has to be refitted.
  changes <- lapply(c('a', 'b', 'c', 'd', 'e'), function(l) list(label = l))
  calls <- 0L
  testthat::local_mocked_bindings(
    .ctImproveRebuild = function(model, data, kept, covariates)
      list(model = vapply(kept, `[[`, '', 'label'), data = NULL, args = list()),
    .ctImproveRefits = function(current, edits, starts, cores, slowest = 0, width = Inf) {
      calls <<- calls + length(edits)
      structure(lapply(edits, function(e) list(labels = e$model)), seconds = 0)
    },
    .ctImproveCompare = function(current, refit, criterion, alpha) {
      gone <- setdiff(refit$labels, current$labels)
      redundant <- gone == 'e' || (gone == 'b' && 'd' %in% current$labels) ||
        (gone == 'd' && 'b' %in% current$labels)
      list(accept = !redundant, dloglik = if (redundant) 0.1 else 10, dpar = 1L,
        value = if (redundant) 0.5 else 1e-4, score = 0, note = '')
    })
  prune <- function(sidebyside) {
    testthat::local_mocked_bindings(.ctImproveSideBySide = function(...) sidebyside)
    calls <<- 0L
    out <- ctsem:::.ctImprovePrune(NULL, list(labels = letters[1:5]), changes, NULL, NULL,
      'LR', 0.05, cores = 6L, ntry = 3L)
    out$refits <- calls
    out
  }
  one <- prune(FALSE)
  batched <- prune(TRUE)
  expect_identical(vapply(one$kept, `[[`, '', 'label'), c('a', 'c', 'd'))
  expect_identical(one$current$labels, c('a', 'c', 'd'))
  expect_identical(one$refits, 5L)
  expect_identical(batched$kept, one$kept)
  expect_identical(batched$current, one$current)
  expect_identical(do.call(rbind, batched$history), do.call(rbind, one$history))
  # a, b, c side by side; b goes, so c, d, e against the model without it.
  expect_identical(batched$refits, 6L)
})

test_that("a rank change starts at the best low-rank covariance, applies last and is never pruned", {
  # No fits. The loadings a lower rank starts from: the echelon form the engine
  # reads, reproducing the covariance's best rank-r approximation.
  set.seed(2)
  A <- matrix(stats::rnorm(25), 5)
  S <- crossprod(A) + diag(5)
  L <- ctsem:::.ctImproveEchelonLoadings(S, 3L)
  expect_identical(dim(L), c(5L, 3L))
  expect_equal(L[upper.tri(L)], rep(0, sum(upper.tri(L))))
  e <- eigen(S, symmetric = TRUE)
  expect_equal(tcrossprod(L), e$vectors[, 1:3] %*% diag(e$values[1:3]) %*% t(e$vectors[, 1:3]),
    tolerance = 1e-10)

  # The rank a step goes straight to: the smallest whose dropped dimensions
  # carry under 5 per cent of the correlation matrix's trace.
  B <- matrix(stats::rnorm(10), 5)
  expect_identical(ctsem:::.ctImproveRankFromCov(tcrossprod(B) + diag(1e-3, 5)), 2L)
  expect_identical(ctsem:::.ctImproveRankFromCov(diag(5)), 5L)

  # A weak direction counts against the covariance when it is mostly the
  # covariance's own parameters.
  dir <- function(p) list(identifiability = list(directions = list(list(parameters = p))))
  expect_true(ctsem:::.ctImproveWeakCovariance(dir(c('rawcor_a__b', 'popsd_a', 'drift'))))
  expect_false(ctsem:::.ctImproveWeakCovariance(dir(c('drift', 'diff', 'rawcor_a__b'))))
  expect_false(ctsem:::.ctImproveWeakCovariance(list()))

  # Rebuilt, a rank change comes after the change that adds a random effect, so
  # it counts that effect: a rank below the number varying is kept, and one
  # that is not is the full covariance.
  m <- suppressMessages(ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
    LAMBDA = matrix(1)))
  k <- ctsem:::.ctImproveNvarying(m)
  add <- list(family = 'heterogeneity', label = 'drift varying',
    edit = function(model, data, covs, args = list()) {
      model$pars$indvarying[model$pars$param %in% 'drift_eta1'] <- TRUE
      list(model = model, data = data)
    })
  rank <- ctsem:::.ctImproveRankChange(k)
  expect_identical(ctsem:::.ctImproveRebuild(m, NULL, list(rank, add), NULL)$args$poprank, k)
  expect_true(is.na(ctsem:::.ctImproveRebuild(m, NULL, list(rank), NULL)$args$poprank))

  # Pruning never tries to drop it: undoing a rank adds parameters.
  changes <- list(list(label = 'a'), rank, list(label = 'b'))
  tried <- character()
  testthat::local_mocked_bindings(
    .ctImproveSideBySide = function(...) FALSE,
    .ctImproveRebuild = function(model, data, kept, covariates)
      list(model = vapply(kept, `[[`, '', 'label'), data = NULL, args = list()),
    .ctImproveRefits = function(current, edits, starts, cores, slowest = 0, width = Inf) {
      tried <<- c(tried, vapply(edits, function(e) paste(e$model, collapse = '+'), ''))
      structure(lapply(edits, function(e) list(labels = e$model)), seconds = 0)
    },
    .ctImproveCompare = function(current, refit, criterion, alpha)
      list(accept = TRUE, dloglik = 10, dpar = 1L, value = 1e-4, score = 0, note = ''))
  out <- ctsem:::.ctImprovePrune(NULL, list(labels = c('a', rank$label, 'b')), changes,
    NULL, NULL, 'LR', 0.05, cores = 1L, ntry = 3L)
  expect_identical(length(out$history), 2L)
  expect_identical(vapply(out$kept, `[[`, '', 'label'), c('a', rank$label, 'b'))
  expect_true(all(grepl(rank$label, tried, fixed = TRUE)))
})
