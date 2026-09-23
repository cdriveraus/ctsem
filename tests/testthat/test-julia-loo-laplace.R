# ctLOO on intoverpop = 'laplace' fits: what K fold scores, and the refit-free
# PSIS method.
#
# The model is a linear Gaussian random intercept with every dynamic parameter
# fixed, so the Laplace approximation is exact and only two population
# parameters are free -- the intercept mean and its population sd. That is
# what makes the oracles below possible: the exact leave-one-subject-out
# density is a two-dimensional integral over the population parameters, and
# the exact leave-one-row-out density a one-dimensional integral over each
# subject's intercept, and both are done here on a grid rather than by the
# importance sampling under test.

.lool_model <- function() {
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
.lool_data <- function(nsubjects = 16, nobs = 5, seed = 11) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(nsubjects), function(i) {
    intercept <- stats::rnorm(1, 1, 1.2)
    state <- stats::rnorm(1, 0, 0.5)
    y <- numeric(nobs)
    for (t in seq_len(nobs)) {
      if (t > 1) {
        state <- exp(-0.5) * state +
          stats::rnorm(1, 0, sqrt(0.36 / 1 * (1 - exp(-1))))
      }
      y[t] <- state + intercept + stats::rnorm(1, 0, 0.4)
    }
    data.frame(id = i, time = seq_len(nobs) - 1, Y1 = y)
  }))
}

.lool_cache <- new.env(parent = emptyenv())
.lool_fit <- function() {
  if (is.null(.lool_cache$fit)) {
    .lool_cache$fit <- suppressMessages(ctFit(.lool_data(), .lool_model(),
      backend = "julia", intoverpop = "laplace", verbose = 0))
  }
  .lool_cache$fit
}

.lool_unit_terms <- function(fit, pars = fit$estimate$raw) {
  ctsem:::.ctBackendLaplaceUnitTerms(ctsem:::.ctBackendSpec(fit),
    matrix(as.numeric(pars), ncol = 1L))
}

test_that("K fold row scoring on a laplace fit keeps the held-out row out of the mode", {
  skip_without_julia()
  fit <- .lool_fit()
  data <- ctsem:::.ctBackendSpec(fit)$data
  set.seed(1); before <- suppressMessages(ctLOO(fit, folds = 4, refit = FALSE,
    subjectwise = FALSE, cores = 1))
  expect_identical(before$scoring,
    "rows at random effect modes from the fold's training rows")

  # Shift one observation a long way. It sits in exactly one fold; in that
  # fold's scoring, the rows of the same subject that come *before* it depend
  # on it only through the random-effect mode, so if the mode was solved from
  # the training rows they cannot move at all.
  target <- 3L
  fold <- which(vapply(before$foldrows, function(x) target %in% x, logical(1)))
  shifted <- data; shifted$Y1[target] <- shifted$Y1[target] + 5
  refit <- ctsem:::.ctFitReplaceData(fit, shifted)
  set.seed(1); after <- suppressMessages(ctLOO(refit, folds = 4, refit = FALSE,
    subjectwise = FALSE, cores = 1))
  expect_equal(unname(after$foldrows), unname(before$foldrows))
  earlier <- which(data$id == data$id[target] & seq_len(nrow(data)) < target)
  expect_equal(after$LogLikRowFolds[[fold]][earlier],
    before$LogLikRowFolds[[fold]][earlier], tolerance = 1e-10)

  # And the full-data scoring, which does see the row, does move: the check
  # above is not passing because nothing depends on anything.
  expect_false(isTRUE(all.equal(as.numeric(after$insampleLogLikRow)[earlier],
    as.numeric(before$insampleLogLikRow)[earlier])))
  expect_lt(before$outsampleLogLik, before$insampleLogLik)
})

test_that("K fold subject folds on a laplace fit score each subject's marginal", {
  skip_without_julia()
  fit <- .lool_fit()
  set.seed(2); out <- suppressMessages(ctLOO(fit, folds = 4, refit = FALSE,
    subjectwise = TRUE, cores = 1))
  marginal <- .lool_unit_terms(fit)$unit_loglik[, 1]

  # At fixed parameters a one-level subject's out-of-sample value is its own
  # Laplace marginal: no mode fitted to its data enters.
  expect_equal(out$outsampleLogLikSubject, marginal, tolerance = 1e-10)
  expect_equal(out$outsampleLogLik, sum(marginal), tolerance = 1e-10)
  expect_true(all(is.na(out$outsampleLogLikRow)))
  expect_identical(out$scoring,
    "subject marginal at the fold's parameters; rows not scored")

  # The quantity that used to be reported -- that subject's rows at a mode
  # fitted to those same rows -- is not the marginal, and is optimistic.
  rowsubject <- ctsem:::.ctFitRowSubject(fit)
  conditional <- vapply(seq_along(marginal), function(i)
    sum(out$insampleLogLikRow[1, rowsubject == i], na.rm = TRUE), numeric(1))
  expect_true(all(conditional > marginal))

  # Changing one subject's data changes its score through the marginal and
  # nothing else, and leaves every other subject's alone.
  data <- ctsem:::.ctBackendSpec(fit)$data
  data$Y1[data$id == 5] <- data$Y1[data$id == 5] + 2
  moved <- ctsem:::.ctFitReplaceData(fit, data)
  set.seed(2); out2 <- suppressMessages(ctLOO(moved, folds = 4, refit = FALSE,
    subjectwise = TRUE, cores = 1))
  expect_equal(out2$outsampleLogLikSubject, .lool_unit_terms(moved)$unit_loglik[, 1],
    tolerance = 1e-10)
  expect_equal(out2$outsampleLogLikSubject[-5], out$outsampleLogLikSubject[-5],
    tolerance = 1e-10)
})

test_that("K fold with refitting runs on a laplace fit at both fold types", {
  skip_without_julia()
  fit <- .lool_fit()
  set.seed(3); rows <- suppressMessages(ctLOO(fit, folds = 2, cores = 1,
    subjectwise = FALSE))
  expect_false(anyNA(rows$outsampleLogLikRow[!is.na(rows$insampleLogLikRow[1, ])]))
  expect_lt(rows$outsampleLogLik, rows$insampleLogLik)
  set.seed(3); subj <- suppressMessages(ctLOO(fit, folds = 2, cores = 1,
    subjectwise = TRUE))
  expect_false(anyNA(subj$outsampleLogLikSubject))
  expect_lt(subj$outsampleLogLik, subj$insampleLogLik)
})

test_that("PSIS leave-one-subject-out matches the exact integral over the population", {
  skip_without_julia()
  skip_if_not_installed("loo")
  fit <- .lool_fit()
  est <- as.numeric(fit$estimate$raw)
  expect_length(est, 2L)

  set.seed(4)
  psis <- suppressWarnings(ctLOO(fit, method = "psis", ndraws = 1000))
  expect_s3_class(psis, "ctLOOpsis")
  expect_equal(nrow(psis$pointwise), 16L)
  expect_equal(psis$ndropped, 0L)

  # Exact: log p(y_i | y_-i) = log int p(y | th) dth - log int p(y_-i | th) dth,
  # on a grid in the covariance's own axes wide enough that the tails beyond it
  # do not matter. The unit terms are exact marginals for this model, and the
  # log prior is the objective's own, so this is the target PSIS estimates.
  eig <- eigen(fit$estimate$cov, symmetric = TRUE)
  axis <- seq(-9, 9, length.out = 91)
  nodes <- as.matrix(expand.grid(a = axis, b = axis))
  grid <- est + eig$vectors %*% (sqrt(eig$values) * t(nodes))
  terms <- ctsem:::.ctBackendLaplaceUnitTerms(ctsem:::.ctBackendSpec(fit), grid)
  expect_true(all(is.finite(terms$value)))
  lse <- ctsem:::.ctLogSumExp
  exact <- vapply(seq_len(16), function(i)
    lse(terms$value) - lse(terms$value - terms$unit_loglik[i, ]), numeric(1))

  # Absolute, elementwise: the largest per-subject error, then the total.
  expect_lt(max(abs(psis$pointwise$elpd_loo - exact)), 0.02)
  expect_lt(abs(psis$elpd_loo - sum(exact)), 0.1)
  expect_true(all(psis$pointwise$pareto_k < 0.7))
  # And it is not the plug-in: integrating over the population is worse than
  # scoring at the estimate, by an amount that is not Monte Carlo noise.
  expect_lt(psis$elpd_loo, sum(.lool_unit_terms(fit)$unit_loglik) - 0.5)
})

test_that("PSIS leave-one-row-out matches a grid over each subject's intercept", {
  skip_without_julia()
  skip_if_not_installed("loo")
  fit <- .lool_fit()
  est <- as.numeric(fit$estimate$raw)
  set.seed(5)
  psis <- suppressWarnings(ctLOO(fit, method = "psis", subjectwise = FALSE,
    ndraws = 1000))
  expect_identical(psis$level, "row")
  expect_equal(nrow(psis$pointwise), nrow(ctsem:::.ctBackendSpec(fit)$data))

  # Every row's one-step-ahead log likelihood with every subject's standardised
  # intercept at each grid value, in one Julia call; the draws and weights under
  # test play no part.
  JuliaConnectoR::juliaEval(paste(
    "function _ctsem_test_loo_grid(obj, theta, grid)",
    "  n = length(obj.objective.subject_objectives)",
    "  out = zeros(0, 0)",
    "  for (g, u) in enumerate(grid)",
    "    sv = ContinuousTimeSEM.ctsem_laplace_subject_values(obj, theta, fill(u, n))",
    "    tr = ContinuousTimeSEM.ctsem_kalman(obj.objective, sv;",
    "      subject_matrices=false, fields=[\"llrow\"])",
    "    g == 1 && (out = zeros(length(tr.llrow), length(grid)))",
    "    out[:, g] .= tr.llrow",
    "  end",
    "  return out",
    "end", sep = "\n"))
  grid <- seq(-9, 9, by = 0.004)
  ll <- JuliaConnectoR::juliaCall("_ctsem_test_loo_grid",
    ctsem:::.ctJuliaObjective(fit), as.numeric(est), grid)
  ll <- matrix(as.numeric(ll), ncol = length(grid))
  rowsubject <- ctsem:::.ctFitRowSubject(fit)
  lse <- ctsem:::.ctLogSumExp
  prior <- -grid^2 / 2
  exact <- vapply(seq_len(nrow(ll)), function(r) {
    total <- prior + colSums(ll[rowsubject == rowsubject[r], , drop = FALSE])
    lse(total) - lse(total - ll[r, ])
  }, numeric(1))

  expect_lt(max(abs(psis$pointwise$elpd_loo - exact[psis$pointwise$row])), 0.01)
  expect_true(all(psis$pointwise$pareto_k < 0.7))
})

test_that("PSIS is refused by name where it does not apply", {
  skip_without_julia()
  fit <- .lool_fit()
  expect_error(ctLOO(fit, method = "psis", folds = 3), "folds is not used")
  expect_error(ctLOO(fit, method = "psis", refit = FALSE, leaveOutN = 2),
    "leaveOutN, refit are not used")
  expect_error(ctLOO(fit, ndraws = 100), "ndraws is used only by method = 'psis'")

  estonly <- fit
  estonly$estimate$cov <- NULL
  expect_error(ctLOO(estonly, method = "psis"), "carries no covariance")

  # A julia fit without Laplace random effects.
  fixed <- .lool_model()
  fixed$pars$indvarying <- FALSE
  plain <- suppressMessages(ctFit(.lool_data(nsubjects = 4), fixed,
    backend = "julia", verbose = 0,
    optimcontrol = list(estonly = TRUE)))
  expect_error(ctLOO(plain, method = "psis"), "intoverpop = 'laplace'")
  # A stan fit.
  expect_error(ctLOO(ctstantestfit, method = "psis"), "intoverpop = 'laplace'")
})
