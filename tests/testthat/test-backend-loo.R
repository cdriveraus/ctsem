# K-fold cross validation for backend='julia' fits.
#
# `ctLOO()` is the one function in this family that could not simply be reached
# through the shared accessors: the Stan path withholds rows with
# `standata$dokalmanrows` and refits with `stanoptimis`, and the julia backend
# has neither. It withholds at the *data* level instead -- a held-out row's
# manifests are set to NA, which is a row the filter propagates through without
# an update and without a likelihood contribution -- and refits with the
# engine's own optimizer.
#
# The tests below are the two that would catch a mistake in that translation:
# an identity that only holds if withholding actually withholds, and a
# comparison against Stan on the same data.

.loo_data <- function() {
  set.seed(3)
  times <- c(0, .6, 1.3, 2.1, 3.0, 3.8)
  drift <- -0.7; diffusion <- 0.5
  do.call(rbind, lapply(seq_len(24), function(i) {
    state <- stats::rnorm(1, 0, .6)
    y <- numeric(length(times))
    for (t in seq_along(times)) {
      if (t > 1) {
        dt <- times[t] - times[t - 1]
        state <- exp(drift * dt) * state +
          stats::rnorm(1, 0, diffusion * sqrt((1 - exp(2 * drift * dt)) / (-2 * drift)))
      }
      y[t] <- state + stats::rnorm(1, 0, .3) + 1.1
    }
    data.frame(id = i, time = times, Y1 = y)
  }))
}

.loo_model <- function() {
  suppressWarnings(ctModel(type = "ct", LAMBDA = matrix(1, 1, 1),
    DRIFT = matrix("drift", 1, 1), DIFFUSION = matrix("diff", 1, 1),
    MANIFESTVAR = matrix("mvar", 1, 1),
    MANIFESTMEANS = matrix("mmean||FALSE", 1, 1),
    T0VAR = matrix("t0v", 1, 1), T0MEANS = matrix(0, 1, 1),
    CINT = matrix(0, 1, 1)))
}

test_that("ctLOO runs on a julia fit and reports the same structure as Stan's", {
  skip_without_julia()
  data <- .loo_data()
  fit <- suppressMessages(ctFit(data, .loo_model(), backend = "julia", verbose = 0,
    optimcontrol = list(estonly = TRUE)))

  set.seed(4)
  out <- suppressMessages(ctLOO(fit, folds = 3, cores = 1))

  expect_true(all(c("foldrows", "foldpars", "insampleLogLikRow", "LogLikRowFolds",
    "outsampleLogLikRow", "insampleLogLik", "outsampleLogLik",
    "insampleRowwiseEntropy", "outsampleRowwiseEntropy") %in% names(out)))
  expect_equal(length(out$foldrows), 3L)
  expect_equal(ncol(out$foldpars), 3L)
  expect_equal(nrow(out$foldpars), length(fit$estimate$raw))
  # Stan returns this as a one-row matrix; so does this.
  expect_equal(dim(out$insampleLogLikRow), c(1L, nrow(data)))
  expect_equal(length(out$outsampleLogLikRow), nrow(data))
  # Every row is held out by exactly one fold, so every row has an out-of-sample
  # likelihood -- a fold that silently dropped rows would leave NAs here.
  expect_false(anyNA(out$outsampleLogLikRow))
  expect_equal(sort(unname(unlist(out$foldrows))), seq_len(nrow(data)))

  # Out of sample is worse than in sample. This is the whole point of the
  # exercise, and it is only true if the held-out rows were genuinely withheld
  # from the refit rather than merely re-scored afterwards.
  expect_lt(out$outsampleLogLik, out$insampleLogLik)
})

test_that("without refitting, the out-of-sample likelihoods are the in-sample ones", {
  skip_without_julia()
  # The sharpest available check on the bookkeeping: with `refit=FALSE` every
  # fold scores the same parameters against the same full data, so assembling
  # the out-of-sample vector fold by fold must reconstruct the in-sample one
  # exactly. Any error in which rows a fold owns shows up here as a mismatch,
  # with no optimizer noise to hide behind.
  data <- .loo_data()
  fit <- suppressMessages(ctFit(data, .loo_model(), backend = "julia", verbose = 0,
    optimcontrol = list(estonly = TRUE)))

  set.seed(5)
  out <- suppressMessages(ctLOO(fit, folds = 4, cores = 1, refit = FALSE))
  insample <- as.numeric(out$insampleLogLikRow)
  outsample <- as.numeric(out$outsampleLogLikRow)
  observed <- !is.na(insample)
  expect_equal(outsample[observed], insample[observed])
})

test_that("julia and Stan cross-validate to the same answer", {
  skip_if_not_installed("rstan")
  skip_without_julia()
  data <- .loo_data()
  model <- .loo_model()

  julia_fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0,
    optimcontrol = list(estonly = TRUE)))
  stan_fit <- suppressMessages(ctFit(data, model, backend = "stan", optimize = TRUE,
    optimcontrol = list(carefulfit = FALSE, stochastic = FALSE), cores = 1, verbose = 0))

  # `refit = FALSE` deliberately: what this test is for is the *translation* --
  # which rows each fold owns, how they are withheld, and how the pieces are
  # assembled back into row and subject summaries. Refitting would add two
  # independent optimizations per fold on top of that, and a fold that withholds
  # a third of the data is exactly where this backend's known weakness bites:
  # neither `ctsem_optimize` nor `stanoptimis` restarts from a flat direction,
  # and on this model they land up to 6 raw units apart in one parameter while
  # differing by ~2% in likelihood. Comparing that would be testing the
  # optimizers' luck. The refit path is covered by the identities above.
  set.seed(6); julia <- suppressMessages(ctLOO(julia_fit, folds = 3, cores = 1,
    refit = FALSE))
  set.seed(6); stan <- suppressMessages(ctLOO(stan_fit, folds = 3, cores = 1,
    refit = FALSE))

  expect_identical(sort(names(julia)), sort(names(stan)))
  # The same folds, because the fold construction is the same expression driven
  # by the same seed.
  expect_equal(unname(julia$foldrows), unname(stan$foldrows))
  expect_equal(as.numeric(julia$insampleLogLikRow),
    as.numeric(stan$insampleLogLikRow), tolerance = 1e-4)
  expect_equal(julia$outsampleLogLik, stan$outsampleLogLik, tolerance = 1e-4)
  expect_equal(julia$outsampleRowwiseEntropy, stan$outsampleRowwiseEntropy,
    tolerance = 1e-4)
  expect_equal(julia$outsampleSubjectwiseLogLikSD,
    stan$outsampleSubjectwiseLogLikSD, tolerance = 1e-4)
})

test_that("subjectwise folds hold out whole subjects", {
  skip_without_julia()
  data <- .loo_data()
  fit <- suppressMessages(ctFit(data, .loo_model(), backend = "julia", verbose = 0,
    optimcontrol = list(estonly = TRUE)))

  set.seed(7)
  out <- suppressMessages(ctLOO(fit, folds = 3, cores = 1, subjectwise = TRUE,
    refit = FALSE))
  rowsubject <- ctsem:::.ctFitRowSubject(fit)
  # No subject may appear in two folds, or the "out of sample" rows would have
  # informed the fit through their own subject's other observations.
  for (fold in out$foldrows) {
    subjects <- unique(rowsubject[fold])
    expect_equal(sort(fold), sort(which(rowsubject %in% subjects)))
  }
  expect_equal(sort(unname(unlist(out$foldrows))), seq_len(nrow(data)))
})
