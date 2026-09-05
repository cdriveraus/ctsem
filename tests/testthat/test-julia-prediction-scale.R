# What `ctKalman()` reports for a non-Gaussian indicator.
#
# The filter's internal quantity is a linear predictor -- a logit, a log rate --
# and reporting it unchanged made every comparison against the data a comparison
# between two different things: a binary model's "probability" reached 4.17
# against data in [0,1], and a count model reported [0.850, 1.346] against
# counts in [0, 27] averaging 5.12.
#
# So the assertions here are about *scale*, and each is written so that the
# linear predictor fails it. Range is the weakest of them and is not enough on
# its own -- a plug-in `inv_logit(ηbar)` is in range too, and is still the wrong
# answer, because an observation is a draw from the marginal distribution with
# the state uncertainty integrated out rather than from the one at the state's
# mean. Where a closed form exists (count) or a one-dimensional integral is
# cheap (binary) the reported value is checked against it directly, and its
# distance from the plug-in is checked too, so that swapping one for the other
# cannot pass.

.pscale_cache <- new.env(parent = emptyenv())

# Fitting is the expensive part and several assertions want the same fit, so
# each model is fitted once per session and the whole Kalman output kept.
.pscale <- function(name) {
  if (!is.null(.pscale_cache[[name]])) return(.pscale_cache[[name]])
  built <- .pscale_case(name)
  fit <- suppressWarnings(suppressMessages(ctFit(built$d, built$m,
    backend = "julia", intoverpop = "augmented",
    optimcontrol = list(estonly = TRUE))))
  out <- ctBackendKalman(fit, subjects = "all")
  .pscale_cache[[name]] <- out
  out
}

# One latent AR process, observed through whichever indicators the case needs.
# MANIFESTMEANS and LAMBDA are fixed at 0 and 1 throughout, so the linear
# predictor is the latent state itself: that is what lets the checks below
# recompute the marginal moments in R without extracting any parameters.
.pscale_latent <- function(nsubjects, nobs, seed) {
  set.seed(seed)
  gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "eta", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-0.3), DIFFUSION = matrix(0.8), MANIFESTVAR = matrix(1e-3),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(0), Tpoints = nobs))
  data.frame(ctGenerate(gen, n.subjects = nsubjects, Tpoints = nobs,
    backend = "r"))
}

.pscale_model <- function(names, types, ncategories = NULL, manifestvar = NULL) {
  n <- length(names)
  args <- list(type = "ct", n.latent = 1, n.manifest = n, manifestNames = names,
    latentNames = "eta1", manifesttype = as.integer(types),
    LAMBDA = matrix(1, n, 1), MANIFESTMEANS = matrix(0, n, 1),
    CINT = matrix(0), T0MEANS = matrix(0),
    MANIFESTVAR = diag(if (is.null(manifestvar)) rep(0, n) else manifestvar,
      nrow = n))
  if (!is.null(ncategories)) args$ncategories <- as.integer(ncategories)
  m <- suppressWarnings(suppressMessages(do.call(ctModel, args)))
  m$pars$indvarying <- FALSE
  m
}

.pscale_case <- function(name) {
  invlog <- function(x) 1 / (1 + exp(-x))
  if (name == "gaussian") {
    d <- .pscale_latent(20, 6, 3)
    names(d)[names(d) == "eta"] <- "g1"
    return(list(d = d, m = .pscale_model("g1", 0L, manifestvar = 0.3)))
  }
  if (name == "binary") {
    d <- .pscale_latent(25, 8, 11)
    for (i in 1:2) d[[paste0("b", i)]] <- stats::rbinom(nrow(d), 1, invlog(d$eta))
    d$eta <- NULL
    return(list(d = d, m = .pscale_model(c("b1", "b2"), c(1L, 1L))))
  }
  if (name == "ordinal") {
    d <- .pscale_latent(25, 8, 7)
    tau <- c(-1.0, 0.4, 1.9)   # four categories
    draw <- function(eta) {
      cum <- sapply(seq_along(tau), function(k) invlog(tau[k] - eta))
      p <- cbind(cum, 1)
      p <- cbind(p[, 1, drop = FALSE], t(apply(p, 1, diff)))
      apply(p, 1, function(pr) sample.int(length(pr), 1, prob = pmax(pr, 0)))
    }
    for (i in 1:2) d[[paste0("o", i)]] <- draw(d$eta)
    d$eta <- NULL
    return(list(d = d,
      m = .pscale_model(c("o1", "o2"), c(2L, 2L), ncategories = c(4L, 4L))))
  }
  if (name == "count") {
    d <- .pscale_latent(25, 8, 5)
    for (i in 1:2) d[[paste0("c", i)]] <- stats::rpois(nrow(d), exp(d$eta))
    d$eta <- NULL
    return(list(d = d, m = .pscale_model(c("c1", "c2"), c(3L, 3L))))
  }
  if (name == "mixed") {
    d <- .pscale_latent(25, 8, 19)
    d$b1 <- stats::rbinom(nrow(d), 1, invlog(d$eta))
    d$g1 <- d$eta + stats::rnorm(nrow(d), 0, 0.5)
    d$eta <- NULL
    # A measurement variance on the Gaussian row only: the binary row's
    # dispersion comes from its link, and ctFit fixes that cell anyway.
    return(list(d = d, m = .pscale_model(c("b1", "g1"), c(1L, 0L),
      manifestvar = c(0, 0.5))))
  }
  stop("unknown case")
}

.pscale_vec <- function(k, element, column) as.numeric(k[[element]][1, , column])
.pscale_cov <- function(k, element, i, j) as.numeric(k[[element]][1, , i, j])
.pscale_obs <- function(k, column) {
  v <- k$y[, column]
  v[!is.na(v)]
}

test_that("a Gaussian model reports exactly what it reported before", {
  skip_without_julia()
  k <- .pscale("gaussian")
  # LAMBDA is 1 and MANIFESTMEANS 0, so the untouched Gaussian path makes the
  # observation estimate the state itself, to the last bit. A link applied to a
  # `manifesttype = 0` row would show up here.
  for (element in c("prior", "upd", "smooth")) {
    expect_lt(max(abs(.pscale_vec(k, paste0("y", element), 1) -
      .pscale_vec(k, paste0("eta", element), 1))), 1e-14)
  }
  # And `ycov = Jy etacov Jy' + MANIFESTcov` still adds the measurement variance
  # on top of the propagated state covariance -- one constant, positive, and the
  # same at every row.
  gap <- .pscale_cov(k, "ypriorcov", 1, 1) - .pscale_cov(k, "etapriorcov", 1, 1)
  expect_gt(mean(gap), 0)
  expect_lt(stats::sd(gap), 1e-12)
})

test_that("a binary model reports the marginal probability", {
  skip_without_julia()
  k <- .pscale("binary")
  observed <- .pscale_obs(k, 1)
  expect_true(all(observed %in% c(0, 1)))

  for (element in c("yprior", "yupd", "ysmooth")) {
    for (column in 1:2) {
      v <- .pscale_vec(k, element, column)
      expect_true(all(is.finite(v)))
      expect_gte(min(v), 0)
      expect_lte(max(v), 1)
    }
  }
  # A probability on the right scale predicts the observed proportion.
  expect_lt(abs(mean(.pscale_vec(k, "yprior", 1)) - mean(observed)), 0.1)

  # The Bernoulli variance of the reported probability, exactly.
  p <- .pscale_vec(k, "yprior", 1)
  expect_lt(max(abs(.pscale_cov(k, "ypriorcov", 1, 1) - p * (1 - p))), 1e-12)

  # Against the integral itself. LAMBDA is 1 and MANIFESTMEANS 0, so the
  # predictor is `N(etaprior, etapriorcov)` and the reported probability has to
  # be `E[inv_logit]` under it.
  mu <- .pscale_vec(k, "etaprior", 1)
  v <- .pscale_cov(k, "etapriorcov", 1, 1)
  marginal <- mapply(function(mui, vi) {
    if (vi < 1e-10) return(1 / (1 + exp(-mui)))
    stats::integrate(function(e) stats::dnorm(e, mui, sqrt(vi)) / (1 + exp(-e)),
      lower = mui - 10 * sqrt(vi), upper = mui + 10 * sqrt(vi))$value
  }, mu, v)
  expect_lt(max(abs(p - marginal)), 1e-6)
  # And far enough from the plug-in that reporting `inv_logit(ηbar)` instead --
  # which is what Stan reports, and is in range, and is the wrong quantity --
  # fails here.
  expect_gt(max(abs(p - 1 / (1 + exp(-mu)))), 0.005)
})

test_that("a count model reports the marginal rate in closed form", {
  skip_without_julia()
  k <- .pscale("count")
  observed <- .pscale_obs(k, 1)
  expect_true(all(observed >= 0))

  for (element in c("yprior", "yupd", "ysmooth")) {
    for (column in 1:2) {
      v <- .pscale_vec(k, element, column)
      expect_true(all(is.finite(v)))
      expect_gt(min(v), 0)
      expect_lt(max(v), 10 * max(observed) + 10)
    }
  }
  expect_lt(abs(mean(.pscale_vec(k, "yprior", 1)) - mean(observed)),
    0.3 * mean(observed) + 0.1)

  # `E[y] = exp(ηbar + s²/2)` exactly, and MANIFESTMEANS is fixed at zero, so
  # `log(y) - etaprior - etapriorcov/2` is zero at every row and on both
  # indicators. The plug-in `exp(ηbar)` leaves `etapriorcov/2` behind, and that
  # term varies from row to row -- which the line after this one checks, so the
  # zero above is not being met by a term that was zero anyway.
  mu <- .pscale_vec(k, "etaprior", 1)
  v <- .pscale_cov(k, "etapriorcov", 1, 1)
  offsets <- c(log(.pscale_vec(k, "yprior", 1)) - mu - v / 2,
    log(.pscale_vec(k, "yprior", 2)) - mu - v / 2)
  expect_lt(max(abs(offsets)), 1e-9)
  expect_gt(stats::sd(v), 1e-6)

  # Overdispersed relative to a Poisson by exactly the rate's own variance.
  rate <- .pscale_vec(k, "yprior", 1)
  expect_lt(max(abs(.pscale_cov(k, "ypriorcov", 1, 1) -
    (rate + rate^2 * expm1(v)))), 1e-9)
})

test_that("an ordinal model reports the expected category", {
  skip_without_julia()
  k <- .pscale("ordinal")
  observed <- .pscale_obs(k, 1)
  K <- 4

  for (element in c("yprior", "yupd", "ysmooth")) {
    for (column in 1:2) {
      v <- .pscale_vec(k, element, column)
      expect_true(all(is.finite(v)))
      expect_gte(min(v), 1)
      expect_lte(max(v), K)
    }
  }
  expect_lt(abs(mean(.pscale_vec(k, "yprior", 1)) - mean(observed)), 0.4)

  # Every distribution on `1..K` with mean `m` has variance at most
  # `(m-1)(K-m)`, attained only by mass at the two ends. A mean and a variance
  # that break that do not come from one distribution over the categories --
  # which is exactly what a Gaussian variance carried over from the linear
  # predictor is.
  mean_cat <- .pscale_vec(k, "yprior", 1)
  variance <- .pscale_cov(k, "ypriorcov", 1, 1)
  expect_gte(min(variance), 0)
  expect_true(all(variance <= (mean_cat - 1) * (K - mean_cat) + 1e-8))
  expect_gt(max(variance), 0.1)   # not degenerate, so the bound is doing work
})

test_that("a non-Gaussian row leaves the Gaussian rows beside it alone", {
  skip_without_julia()
  k <- .pscale("mixed")
  # Column 1 is binary, column 2 Gaussian.
  binary <- .pscale_vec(k, "yprior", 1)
  expect_gte(min(binary), 0)
  expect_lte(max(binary), 1)

  # The Gaussian row is untouched: still the state, still the state covariance
  # plus a constant measurement variance.
  expect_lt(max(abs(.pscale_vec(k, "yprior", 2) -
    .pscale_vec(k, "etaprior", 1))), 1e-14)
  gap <- .pscale_cov(k, "ypriorcov", 2, 2) - .pscale_cov(k, "etapriorcov", 1, 1)
  expect_gt(mean(gap), 0)
  expect_lt(stats::sd(gap), 1e-12)

  # The cross-covariance between them is dropped rather than carried over: one
  # is a probability and the other a manifest value, and `Jy P Jy'` is a
  # covariance between linear predictors, so the pair have no shared scale.
  for (element in c("ypriorcov", "yupdcov", "ysmoothcov")) {
    expect_true(all(.pscale_cov(k, element, 1, 2) == 0))
    expect_true(all(.pscale_cov(k, element, 2, 1) == 0))
  }
})

test_that("indicators of the same non-Gaussian kind share no covariance either", {
  skip_without_julia()
  for (name in c("binary", "ordinal", "count")) {
    k <- .pscale(name)
    for (element in c("ypriorcov", "yupdcov", "ysmoothcov")) {
      expect_true(all(.pscale_cov(k, element, 1, 2) == 0))
      expect_true(all(.pscale_cov(k, element, 2, 1) == 0))
    }
  }
})
