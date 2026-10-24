# The coreset importance sampling evaluates on the Laplace route (R/ctCoreset.R,
# laplace_coreset.jl, coreset_sampler.jl): the engine's subset evaluation
# against the full one, and the sampler's density with a coreset behind it
# against the full density. The choice on sums it can compress is the
# engine's own test (test_coreset_sampler.jl). One fit, the shared Laplace
# fixture at 100 subjects: its seven parameters all move under 'is', so a
# coreset needs about twice that many units, a small enough share of 100 to
# pay.

.coreset_fixture <- function() laplace_fixture(nsub = 100L, tp = 6L, seed = 31L)

test_that("the engine's coreset evaluation is the weighted sum of the units it names", {
  skip_on_cran()
  skip_if_not(isTRUE(ctsem::ctJuliaStatus()$available))
  fit <- .coreset_fixture()
  module <- ctsem:::.ctJuliaModule(fit$model_spec$project)
  obj <- ctsem:::.ctJuliaObjective(fit, hybrid = TRUE)
  est <- as.numeric(fit$estimate$raw)
  set.seed(3)
  X <- sweep(matrix(rnorm(3 * length(est), 0, 0.05), 3), 2, est, "+")
  terms <- ctsem:::.ctBackendJuliaValue(module$ctsem_coreset_terms(obj, t(X)))
  L <- matrix(as.numeric(terms$unit_loglik), ncol = 3)
  N <- nrow(L)
  expect_equal(N, 100L)
  plain <- as.numeric(ctsem:::.ctBackendJuliaValue(module$ctsem_evaluate_batch(obj, t(X))))
  expect_equal(as.numeric(terms$value), plain, tolerance = 1e-10)
  expect_equal(colSums(L), as.numeric(terms$total), tolerance = 1e-10)
  units <- function(x) ctsem:::.ctJuliaVector(as.integer(x))
  weights <- function(x) ctsem:::.ctJuliaVector(as.numeric(x))
  # Every unit at weight one is the full evaluation, value and gradient.
  expect_equal(as.numeric(ctsem:::.ctBackendJuliaValue(module$ctsem_coreset_batch(obj,
    t(X), units(seq_len(N)), weights(rep(1, N))))), plain, tolerance = 1e-10)
  gfull <- matrix(as.numeric(ctsem:::.ctBackendJuliaValue(
    module$ctsem_evaluate_batch_gradient(obj, t(X)))), ncol = 3)
  gall <- matrix(as.numeric(ctsem:::.ctBackendJuliaValue(
    module$ctsem_coreset_batch_gradient(obj, t(X), units(seq_len(N)), weights(rep(1, N))))),
    ncol = 3)
  expect_equal(gall, gfull, tolerance = 1e-9)
  # A subset: the prior plus its weighted terms, and the weighted scores plus
  # the prior's gradient.
  S <- c(2L, 9L, 23L, 31L)
  w <- c(0.5, 3, 1.25, 7)
  prior <- as.numeric(terms$value) - as.numeric(terms$total)
  expect_equal(as.numeric(ctsem:::.ctBackendJuliaValue(module$ctsem_coreset_batch(obj,
    t(X), units(S), weights(w)))), prior + colSums(L[S, ] * w), tolerance = 1e-10)
  g <- matrix(as.numeric(ctsem:::.ctBackendJuliaValue(
    module$ctsem_coreset_batch_gradient(obj, t(X), units(S), weights(w)))), ncol = 3)
  scores <- as.matrix(ctsem:::.ctJuliaGet(module$ctsem_subject_gradients(obj,
    ctsem:::.ctJuliaVector(X[1, ]), prior = FALSE))$scores)
  expect_equal(g[-1, 1], gfull[-1, 1] - colSums(scores) + colSums(scores[S, ] * w),
    tolerance = 1e-8)
  # One unit crosses the bridge as a vector too, and a repeated unit is refused.
  expect_equal(as.numeric(ctsem:::.ctBackendJuliaValue(module$ctsem_coreset_batch(obj,
    t(X), units(5L), weights(2)))), prior + 2 * L[5, ], tolerance = 1e-10)
  expect_error(module$ctsem_coreset_batch(obj, t(X), units(c(1L, 1L)), weights(c(1, 1))))
  # A gradient batch that teaches: the terms, the batch gradient, and each
  # unit's score at the columns asked for, which sum to the gradient of the
  # units' total.
  tg <- ctsem:::.ctBackendJuliaValue(module$ctsem_coreset_terms_gradient(obj, t(X),
    ctsem:::.ctJuliaVector(c(1L, 3L))))
  expect_equal(as.numeric(tg$value), plain, tolerance = 1e-10)
  expect_equal(as.numeric(tg$total), as.numeric(terms$total), tolerance = 1e-10)
  expect_equal(matrix(as.numeric(tg$unit_loglik), ncol = 3), L, tolerance = 1e-10)
  expect_equal(matrix(as.numeric(tg$gradient), ncol = 3), gfull[-1, ], tolerance = 1e-9)
  sc <- array(as.numeric(tg$scores), c(N, length(est), 2))
  expect_equal(sc[, , 1], scores, tolerance = 1e-8)
  scores3 <- as.matrix(ctsem:::.ctJuliaGet(module$ctsem_subject_gradients(obj,
    ctsem:::.ctJuliaVector(X[3, ]), prior = FALSE))$scores)
  expect_equal(sc[, , 2], scores3, tolerance = 1e-8)
})

test_that("the sampler's density with a coreset behind it follows the full density", {
  skip_on_cran()
  skip_if_not(isTRUE(ctsem::ctJuliaStatus()$available))
  fit <- .coreset_fixture()
  full <- ctsem:::.ctBackendLpgFunc(fit, gradient = FALSE)
  lpg <- ctsem:::.ctCoresetLpg(fit, full)
  expect_false(identical(lpg, full))
  expect_identical(ctsem:::.ctCoresetLpg(fit, full, list(coreset = FALSE)), full)
  # A density on another objective, batch attributes and all, is never
  # replaced by the fit's unit terms.
  other <- full
  attr(other, "unitsum") <- NULL
  expect_identical(ctsem:::.ctCoresetLpg(fit, other), other)
  est <- as.numeric(fit$estimate$raw)
  cov <- as.matrix(fit$estimate$cov)
  set.seed(9)
  draw <- function(n) ctsem:::ctOptimNormalDraws(est, cov * 2, n)
  # A full gradient batch teaches a coreset its draws' scores, and a value
  # batch at the same draws -- DAIS asks for both -- costs nothing more and
  # is the full density's.
  taught <- ctsem:::.ctCoresetLpg(fit, full)
  xg <- draw(8)
  gg <- attr(taught, "gradbatch")(xg)
  gt <- attr(full, "gradbatch")(xg)
  expect_equal(gg$value, gt$value, tolerance = 1e-10)
  expect_equal(gg$gradient, gt$gradient, tolerance = 1e-9)
  expect_equal(attr(taught, "cost")(), 16)
  expect_equal(attr(taught, "batch")(xg), attr(full, "batch")(xg), tolerance = 1e-10)
  expect_equal(attr(taught, "cost")(), 16)
  rec0 <- attr(taught, "coreset")()
  expect_equal(rec0$score_draws, 4L)
  expect_equal(rec0$reused_values, 8L)
  # Rounds with no coreset are the full density, half a round at a time; once
  # one is chosen, the rest of the round runs on it, checked (at most three).
  first <- NULL
  v1 <- NULL
  t1 <- NULL
  for (round in 1:3) {
    x <- draw(100)
    had <- isTRUE(attr(lpg, "coreset")()$used)
    v <- attr(lpg, "batch")(x)
    truth1 <- attr(full, "batch")(x)
    if (!had) expect_equal(v[1:50], truth1[1:50], tolerance = 1e-8)
    near <- truth1 > max(truth1) - 15
    expect_lt(max(abs((v - truth1)[near])), 1)
    first <- rbind(first, x)
    v1 <- c(v1, v)
    t1 <- c(t1, truth1)
    if (isTRUE(attr(lpg, "coreset")()$used)) break
  }
  rec <- attr(lpg, "coreset")()
  expect_true(rec$used)
  expect_equal(rec$units, 100L)
  expect_lt(rec$expected_round_cost, 100 + 2 * 30)
  expect_lte(attr(lpg, "cost")(), nrow(first))
  spent <- attr(lpg, "cost")()
  # Rounds on the coreset. A round whose check fails costs a full one (one
  # first coreset in five, over five seeds on this model, against about 20
  # for a round that passes); together they cost well under full rounds.
  second <- NULL
  v2 <- NULL
  truth <- NULL
  for (round in 1:4) {
    x <- draw(100)
    second <- rbind(second, x)
    v2 <- c(v2, attr(lpg, "batch")(x))
    truth <- c(truth, attr(full, "batch")(x))
  }
  # Away from where nothing has weight, the coreset is within a nat.
  near <- truth > max(truth) - 15
  expect_lt(max(abs((v2 - truth)[near])), 1)
  expect_lt(attr(lpg, "cost")() - spent, 300)
  g <- attr(lpg, "gradbatch")(second[1:3, , drop = FALSE])
  gt <- attr(full, "gradbatch")(second[1:3, , drop = FALSE])
  expect_lt(max(abs(g$gradient - gt$gradient)) / max(abs(gt$gradient)), 0.1)
  # The end check: given the coreset's own values it either passes them or
  # returns weights corrected towards the full density's.
  allx <- rbind(first, second)
  logp <- c(v1, v2)
  logw <- logp - max(logp)
  fixed <- attr(lpg, "finalize")(allx, logp, logw)
  rec <- attr(lpg, "coreset")()
  expect_false(is.null(rec$final))
  if (!is.null(fixed)) {
    top <- order(-logw)[1:20]
    tw <- c(t1, truth)[top]
    expect_lt(max(abs((fixed[top] - logw[top]) - (tw - logp[top]))), 1e-6)
  }
})

test_that("importance sampling with a coreset reaches the full density's posterior", {
  skip_on_cran()
  skip_if_not(isTRUE(ctsem::ctJuliaStatus()$available))
  fit <- .coreset_fixture()
  # Each at an effective size of 500. At the default 200 one run's means lie
  # up to 0.22 posterior sds from a long reference (ESS 4,000) on this model
  # and two runs' differ by up to 0.43, with either sampler (seeds 21-24), so
  # the bar below would test the seed; at 500 they differ by 0.08-0.24
  # (seeds 21-26).
  run <- function(control) {
    set.seed(21)
    suppressWarnings(ctsem::ctFitUncertainty(fit, "is", cores = 1,
      control = c(control, list(isESS = 500))))
  }
  a <- run(list(coreset = FALSE))
  b <- run(list())
  ia <- a$uncertainty$details$importance_sampling
  ib <- b$uncertainty$details$importance_sampling
  expect_null(ia$coreset)
  expect_true(isTRUE(ib$coreset$used))
  expect_lt(ib$cost, ia$cost)
  # Within Monte Carlo error of each other, read from the weighted samples.
  expect_equal(b$estimate$se, a$estimate$se, tolerance = 0.25)
  sa <- sqrt(diag(a$uncertainty$imis$covariance))
  expect_lt(max(abs(b$uncertainty$imis$mean - a$uncertainty$imis$mean) / sa), 0.4)
})
