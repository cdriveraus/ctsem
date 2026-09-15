# Count manifest variables: `manifesttype = 3`, Poisson with a log link.
#
# Unlike test-julia-binary.R and test-julia-ordinal.R, there is no closed-form
# single-observation check here, at this level or in the engine's own suite --
# no file under inst/julia/ContinuousTimeSEM/test mentions count, censored, or
# the quadrature functions that carry them (`_binary_moments`, `_binary_mode`,
# `_ekf_binary_update!`) at all. What this file does check: that the model
# accepts and describes the type, that the data checks fire, that the backend
# that cannot fit one refuses instead of returning a number, that the
# reverse-mode gradient is right on a real model, and that a fit recovers
# generating parameters.

.count_data <- function(nsubjects = 40, nobs = 8, nind = 2, seed = 5) {
  set.seed(seed)
  d <- do.call(rbind, lapply(seq_len(nsubjects), function(i) {
    gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
      manifestNames = "eta", latentNames = "eta1", LAMBDA = matrix(1),
      DRIFT = matrix(-0.4), DIFFUSION = matrix(0.6),
      MANIFESTVAR = matrix(1e-6), T0VAR = matrix(1), T0MEANS = matrix(0),
      CINT = matrix(0), MANIFESTMEANS = matrix(0), Tpoints = nobs))
    one <- data.frame(ctGenerate(gen, n.subjects = 1, Tpoints = nobs,
      backend = "r"))
    one$id <- i
    one
  }))
  for (j in seq_len(nind)) {
    d[[paste0("c", j)]] <- stats::rpois(nrow(d), exp(1.2 + d$eta))
  }
  d$eta <- NULL
  d
}

.count_model <- function(nind = 2) {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = nind, manifestNames = paste0("c", seq_len(nind)),
    latentNames = "eta1", manifesttype = rep(3L, nind),
    LAMBDA = matrix(1, nind, 1),
    MANIFESTMEANS = matrix("mm", nind, 1), CINT = matrix(0),
    T0MEANS = matrix(0), MANIFESTVAR = diag(0, nind))))
  m$pars$indvarying <- FALSE
  m
}

# The same model with something for the Laplace route to integrate over.
# `intoverpop='laplace'` refuses a model with no varying parameter, correctly,
# so the two-route comparisons need one.
.count_model_varying <- function(nind = 2) {
  m <- .count_model(nind)
  m$pars$indvarying[m$pars$param %in% "mm"] <- TRUE
  m
}

test_that("ctModel takes manifesttype 3 and needs no categories for it", {
  m <- .count_model()
  expect_equal(unname(m$manifesttype), c(3L, 3L))
  # `ncategories` describes ordinal variables only, and a count must not be
  # made to invent one.
  expect_true(all(m$ncategories == 0L))
  expect_null(m$THRESHOLDS)
  expect_error(suppressWarnings(suppressMessages(ctModel(type = "ct",
    n.latent = 1, n.manifest = 1, manifestNames = "y", latentNames = "eta1",
    manifesttype = 5L, LAMBDA = matrix(1)))), "manifesttype must be")
})

test_that("print names counts and the link they use", {
  described <- paste(utils::capture.output(print(.count_model())),
    collapse = " ")
  expect_match(described, "c1, c2 \\(count, Poisson log link\\)")
})

test_that("count data is checked against the model", {
  d <- .count_data(nsubjects = 5, nobs = 4)
  m <- .count_model()

  negative <- d
  negative$c1[3] <- -1
  expect_error(suppressMessages(ctFit(negative, m, backend = "julia",
    fit = FALSE)), "negative values")

  fractional <- d
  fractional$c1[3] <- 1.5
  expect_error(suppressMessages(ctFit(fractional, m, backend = "julia",
    fit = FALSE)), "whole number")

  # Legal, and worth saying: a constant count identifies no rate.
  constant <- d
  constant$c1 <- 0L
  expect_warning(suppressMessages(ctFit(constant, m, backend = "julia",
    fit = FALSE)), "not identified")
})

test_that("stan refuses counts rather than treating them as Gaussian", {
  d <- .count_data(nsubjects = 5, nobs = 4)
  expect_error(suppressMessages(ctFit(d, .count_model(), backend = "stan",
    optimcontrol = list(estonly = TRUE))), "need backend=\"julia\"")
})

test_that("the adjoint matches forward mode on a count model", {
  skip_without_julia()
  d <- .count_data()
  fit <- suppressWarnings(suppressMessages(ctFit(d, .count_model(),
    backend = "julia", intoverpop = "augmented",
    optimcontrol = list(estonly = TRUE))))
  handle <- structure(fit$model_spec,
    class = c("ctJuliaModel", "ctFitModel"))
  est <- as.numeric(fit$estimate$raw)

  # Anchored at the estimate rather than drawn from zero, which is what the
  # binary and ordinal suites do. Their likelihoods are bounded, so any raw
  # draw gives a representable objective; a count's linear predictor is
  # *exponentiated*, and the same draw from zero puts the Poisson rate near
  # 1e29. Comparing two gradients there compares two numbers that have lost
  # every significant digit -- measured, it reported a 50% disagreement where
  # the same code agrees to 1e-16 anywhere the objective is representable.
  set.seed(11)
  for (trial in 1:3) {
    at <- est + stats::rnorm(length(est), 0, 0.3)
    adjoint <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
      gradient_method = "adjoint")$gradient)
    forward <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
      gradient_method = "forward")$gradient)
    expect_equal(adjoint, forward, tolerance = 1e-8)
  }
})

test_that("a count model recovers what generated it, and laplace does it better", {
  skip_without_julia()
  d <- .count_data()
  m <- .count_model_varying()

  # A random start, deliberately: `inits = NULL` draws `rnorm(npar, 0, .01)`
  # from wherever the session's RNG has reached, so this asks whether the fit
  # depends on where it began. It used to. Over twenty such starts the laplace
  # route reached the optimum 14 times and the augmented route 9, and eleven of
  # the augmented failures reported `converged = TRUE` at a largest gradient of
  # 2.7e4; both are 20/20 now. The two defects behind that are a count
  # intercept carrying a location parameter's `meanscale` (see
  # `.ctModelDefaultFreePar`) and an inner mode solve that reported failure for
  # reaching the limit of double precision (see `_laplace_newton_unit_mode`).
  #
  # So this test is not pinned, and should not be: a pinned start would pass
  # over either defect returning.
  fits <- lapply(c(augmented = "augmented", laplace = "laplace"), function(route) {
    suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
      intoverpop = route, optimcontrol = list(estonly = TRUE))))
  })
  for (route in names(fits)) expect_true(isTRUE(fits[[route]]$optim$converged))

  got <- vapply(fits, function(f) {
    means <- summary(f)$popmeans
    c(mm = unname(means["mm", "mean"]),
      drift = unname(means["drift_eta1", "mean"]))
  }, numeric(2))

  # This test used to assert that the laplace route was the more accurate of
  # the two here, and record the gap:
  #
  #   route      mm (truth 1.2)      drift (truth -0.4)
  #   laplace    1.242   3.5% off    -0.402   0.5% off
  #   augmented  1.716    43% off    -0.219    45% off
  #
  # attributing it to the augmented route carrying the random intercept as a
  # state through a linearised filter. Almost all of it was the count mode
  # solve instead. `_binary_mode` stopped six Newton steps short, which
  # mis-centred the quadrature by an amount that grows with the linear
  # predictor's variance -- and the augmented route's extra state is exactly an
  # extra contribution to that variance, so it was the route being punished.
  # With the mode solved, both converged:
  #
  #   route      mm (truth 1.2)      drift (truth -0.4)
  #   laplace    1.198  0.17% off    -0.357   10.7% off
  #   augmented  1.198  0.17% off    -0.357   10.7% off
  #
  # The remaining drift miss is this fixture's, not a route's: 40 subjects and
  # 8 occasions, and both routes now find the same maximum. Note that laplace's
  # old drift of -0.402 was luck rather than accuracy -- it moved too, and away
  # from the generating value, because the likelihood it maximises changed.
  #
  # So the claim worth defending is no longer an ordering but an agreement, and
  # it is asserted as one. A reappearance of the old gap fails on the second
  # block below, not on a tolerance wide enough to hide it.
  expect_equal(unname(got["mm", "laplace"]), 1.2, tolerance = 0.05)
  expect_equal(unname(got["drift", "laplace"]), -0.4, tolerance = 0.15)
  expect_equal(unname(got["mm", "augmented"]), 1.2, tolerance = 0.05)
  expect_equal(unname(got["drift", "augmented"]), -0.4, tolerance = 0.15)

  # The two routes integrate the same random intercept two different ways, so
  # on a model this size they should reach the same estimate. Elementwise, so a
  # failure names the parameter that moved.
  expect_equal(unname(got[, "laplace"]), unname(got[, "augmented"]),
    tolerance = 1e-3)
})

test_that("the two routes agree on a count model", {
  skip_without_julia()
  d <- .count_data()
  m <- .count_model_varying()
  # A COMMON, FIXED starting point, because the difference below is read as a
  # methodological gap and that reading needs both routes to have started from
  # the same place. `inits = NULL` starts each at rnorm(npar, 0, .01) from
  # whatever RNG state it inherits, so the laplace fit's start depended on how
  # much RNG the augmented fit above it had consumed -- and on this model that
  # matters far more than 0.01 suggests. From fixed starts:
  #
  #   route      zeros          rnorm sd .01     rnorm sd .3
  #   augmented  -1483.480158   -1483.480158     -1483.480158
  #   laplace    -1471.584356   -1471.584356     -1.47e14
  #
  # The augmented route is flat over all three; the laplace route reaches the
  # intended optimum from a small start, runs away from a wider one, and from
  # some sd-0.01 starts lands near -2175. That is what made this fail
  # intermittently, and it is a property of the model rather than of either
  # construction -- both recorded values below are reproduced exactly from
  # zeros, so zeros is the neutral choice and not the one that passes.
  fit <- function(route) {
    spec <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
      intoverpop = route, fit = FALSE)))
    suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
      intoverpop = route, inits = rep(0, ctsem:::.ctBackendNpar(spec)),
      optimcontrol = list(estonly = TRUE))))
  }
  augmented <- fit("augmented")
  laplace <- fit("laplace")
  # The two routes integrate the same random effect differently, so the gap is
  # a real methodological difference and not noise: -1471.584 (laplace) against
  # The two were 11.896 apart, 0.0080 relative, when the count mode solve
  # stopped short; the gap was that defect rather than a difference between the
  # routes, and the augmented route carried more of it because its extra state
  # adds to the linear predictor's variance. Measured now: -1466.970967283
  # (augmented) against -1466.970967298 (laplace), 1.5e-08 apart and 1.0e-11
  # relative.
  #
  # 1e-6 rather than anything tighter because two different integration routes
  # reaching a maximum by different paths have no reason to agree to the last
  # bit, and it is still four orders of magnitude below the gap this is here to
  # catch.
  expect_equal(as.numeric(laplace$estimate$loglik),
    as.numeric(augmented$estimate$loglik), tolerance = 1e-6)
  # No ordering is asserted. It used to be -- laplace the better fit, as the
  # methodology predicts for an integrated random effect against a linearised
  # one -- and with the mode solved the two agree to 1.5e-08, which is optimiser
  # noise and falls either way between runs. Asserting a direction across a gap
  # that small tests the optimiser's last digit, not the methodology.
})

test_that("ctGenerate draws counts rather than continuous values", {
  skip_without_julia()
  gen <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 2, manifestNames = c("c1", "c2"), latentNames = "eta1",
    manifesttype = c(3L, 3L), LAMBDA = matrix(1, 2, 1),
    DRIFT = matrix(-0.4), DIFFUSION = matrix(0.6),
    MANIFESTVAR = diag(0, 2), T0VAR = matrix(1), T0MEANS = matrix(0),
    CINT = matrix(0), MANIFESTMEANS = matrix(1.0, 2, 1), Tpoints = 6)))
  set.seed(4)
  d <- data.frame(ctGenerate(gen, n.subjects = 15, Tpoints = 6,
    backend = "julia"))
  values <- c(d$c1, d$c2)
  values <- values[!is.na(values)]
  expect_true(length(values) > 0)
  # Whole numbers, none negative: the R generator has no measurement link and
  # would return continuous values here, which is the failure this guards.
  expect_true(all(values >= 0))
  expect_true(all(abs(values - round(values)) < 1e-8))
  expect_gt(length(unique(values)), 1)
})

# ---------------------------------------------------------------------------
# A count's dispersion: MANIFESTVAR's diagonal is a log-scale standard
# deviation, so `y | x` is Poisson-lognormal rather than Poisson.

# One count indicator whose linear predictor is a constant plus its dispersion:
# T0VAR and DIFFUSION are zero, so the state contributes nothing and every row
# is an independent draw from the same Poisson-lognormal. That makes the whole
# model's likelihood something this file can compute from the definition.
.count_iid_model <- function(mu = "mu", v = "v") {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 1, manifestNames = "y", latentNames = "eta1",
    manifesttype = 3L, LAMBDA = matrix(1), DRIFT = matrix(-0.5),
    DIFFUSION = matrix(0), T0VAR = matrix(0), T0MEANS = matrix(0),
    CINT = matrix(0), MANIFESTMEANS = matrix(mu), MANIFESTVAR = matrix(v),
    Tpoints = 1)))
  m$pars$indvarying <- FALSE
  m
}

.count_iid_data <- function(n = 300, mu = 1.1, sigma = 0.7, seed = 99) {
  set.seed(seed)
  data.frame(id = seq_len(n), time = 0,
    y = stats::rpois(n, exp(stats::rnorm(n, mu, sigma))))
}

# log P(y) for one observation under a Poisson-lognormal, by a fixed fine grid
# over the linear predictor with log-sum-exp. Nothing here is shared with the
# engine, which is the point: a comparison against the engine's own machinery
# would agree with an inherited mistake.
#
# A grid rather than `integrate()`, which is the obvious choice and is not safe
# over the range this has to cover. Checked against a 2e6-point version of this
# same grid, `integrate()` over `mu +- 40 sigma` returned `-Inf` for `y = 400`
# at `mu = 0, sigma = 2`, was 37 out at `mu = 2, sigma = 0.5`, and was 0.28 to
# 0.46 out for `y` of 0, 1 and 5 at `mu = 6, sigma = 0.5` -- it misses a narrow
# peak inside a wide interval, and gives no sign that it has. The window below
# brackets both the prior and the likelihood's own mode at `log(y + 1/2)`, so
# the peak is always inside it wherever it sits, and 50001 points agree with
# 2000001 to 6.5e-11 over a grid of `mu` in (0, 2, 6), `sigma` in (0.5, 2, 5)
# and `y` in (0, 1, 5, 50, 400, 5000).
.pln_loglik <- function(y, mu, sigma, npoints = 50001) {
  tab <- sort(unique(y))
  lp <- vapply(tab, function(k) {
    lo <- min(mu - 12 * sigma, log(k + 0.5) - 12)
    hi <- max(mu + 12 * sigma, log(k + 0.5) + 12)
    e <- seq(lo, hi, length.out = npoints)
    lf <- k * e - exp(e) - lgamma(k + 1) - (e - mu)^2 / (2 * sigma^2) -
      log(sigma * sqrt(2 * pi))
    mx <- max(lf)
    mx + log(sum(exp(lf - mx)) * (e[2] - e[1]))
  }, numeric(1))
  sum(lp[match(y, tab)])
}

test_that("a count keeps its measurement variance free, as its dispersion", {
  # Binary and ordinal have theirs fixed to a deterministic value, because a
  # normal term on the linear predictor is absorbed into their loadings and
  # thresholds and is not identified. A count has no such freedom -- the
  # Poisson's variance is locked to its mean -- so the parameter is identified
  # and is left alone, with the same default as a Gaussian indicator's.
  d <- .count_iid_data(n = 30)
  free <- suppressWarnings(suppressMessages(ctFit(d, .count_iid_model(),
    backend = "julia", fit = FALSE)))
  fixed <- suppressWarnings(suppressMessages(ctFit(d,
    .count_iid_model(v = 0), backend = "julia", fit = FALSE)))
  expect_equal(max(free$parameter_table$parnumber, na.rm = TRUE),
    max(fixed$parameter_table$parnumber, na.rm = TRUE) + 1)

  # And a binary indicator in the same position does not keep one, so this is
  # a count-specific decision rather than the gate having been removed.
  b <- .count_iid_model()
  b$manifesttype <- 1L
  db <- d; db$y <- as.numeric(db$y > 3)
  bfit <- suppressWarnings(suppressMessages(ctFit(db, b, backend = "julia",
    fit = FALSE)))
  expect_equal(max(bfit$parameter_table$parnumber, na.rm = TRUE),
    max(fixed$parameter_table$parnumber, na.rm = TRUE))
})

test_that("a count's off-diagonal measurement covariance is fixed to zero", {
  # Counts are applied one row at a time, conditionally independent given the
  # state, so the engine reads only the diagonal. A free off-diagonal would be
  # accepted here and ignored there.
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 2, manifestNames = c("c1", "c2"), latentNames = "eta1",
    manifesttype = c(3L, 3L), LAMBDA = matrix(1, 2, 1), CINT = matrix(0),
    T0MEANS = matrix(0), MANIFESTVAR = "auto", Tpoints = 4)))
  m$pars$indvarying <- FALSE
  d <- .count_data(nsubjects = 8, nobs = 4)
  f <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
    fit = FALSE)))
  pt <- f$parameter_table
  off <- pt[pt$matrix %in% "MANIFESTVAR" & pt$row != pt$col, ]
  expect_true(nrow(off) > 0)
  expect_true(all(is.na(off$parnumber)))
  expect_equal(unname(as.numeric(off$value)), rep(0, nrow(off)))
})

test_that("a count's likelihood is the Poisson-lognormal one", {
  skip_without_julia()
  # The test this file said it did not have. Against a reference computed from
  # the definition rather than against the engine's own quadrature, and over a
  # range of predictor sd, because the error that prompted this was invisible
  # below 0.5 and 1766 log units at 1.2.
  d <- .count_iid_data()
  h <- suppressWarnings(suppressMessages(ctFit(d, .count_iid_model(),
    backend = "julia", fit = FALSE)))
  pt <- h$parameter_table
  pt <- pt[!is.na(pt$parnumber), ]
  # The raw vector for a wanted (mu, sigma), by inverting each parameter's own
  # transform as the model states it -- read rather than hardcoded, so this
  # does not quietly test the wrong point if a default transform changes.
  log1p_exp <- function(x) ifelse(x > 30, x, log1p(exp(x)))
  rawfor <- function(target) {
    vapply(seq_len(nrow(pt)), function(i) {
      tfi <- pt$transform[i]
      want <- target[[pt$param[i]]]
      stats::uniroot(function(r) {
        param <- rep(0, nrow(pt))
        param[pt$parnumber[i]] <- r
        eval(parse(text = tfi)) - want
      }, c(-30, 30), extendInt = "yes", tol = 1e-12)$root
    }, numeric(1))[order(pt$parnumber)]
  }
  # The last two points are the corner the mode solve failed in longest: the
  # data are counts around `exp(1.1)`, so evaluating at `mu = 6` makes every
  # observation a heavily over-predicted small one, which is where a start at
  # `log y` rather than `log(y + 1/2)` left the mode as much as 4.8 out. The
  # large sigma is the other half of it -- the error grew with the predictor's
  # variance and was invisible below 0.5.
  for (at in list(c(1.1, 0.3), c(1.1, 0.7), c(1.1, 1.2), c(6, 0.5), c(6, 2))) {
    engine <- ctJuliaEvaluate(h, rawfor(list(mu = at[1], v = at[2])))$value
    expect_equal(engine, .pln_loglik(d$y, at[1], at[2]), tolerance = 1e-5,
      info = paste("mu", at[1], "predictor sd", at[2]))
  }
})

test_that("the adjoint is right with a count dispersion present", {
  skip_without_julia()
  # Two latents, because a one-latent gradient test hid a symmetry bug in this
  # same adjoint for months. The dispersion's cotangent is the new path: it
  # reaches MANIFESTVAR through the predictor's variance rather than through
  # the density, so nothing in the threshold machinery carries it.
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 2,
    n.manifest = 2, manifestNames = c("c1", "c2"),
    latentNames = c("eta1", "eta2"), manifesttype = c(3L, 3L),
    LAMBDA = matrix(c(1, 0, 0, 1), 2, 2), MANIFESTVAR = "diag",
    CINT = matrix(0, 2, 1), MANIFESTMEANS = matrix(c("m1", "m2"), 2, 1),
    Tpoints = 6)))
  m$pars$indvarying <- FALSE
  set.seed(11)
  n <- 20; tp <- 6
  d <- data.frame(id = rep(seq_len(n), each = tp),
    time = rep(seq_len(tp) - 1, times = n),
    c1 = stats::rpois(n * tp, 4), c2 = stats::rpois(n * tp, 3))
  h <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
    intoverpop = "augmented", fit = FALSE)))
  npar <- max(h$parameter_table$parnumber, na.rm = TRUE)
  set.seed(3)
  for (trial in 1:3) {
    at <- stats::rnorm(npar, 0, 0.25)
    adjoint <- as.numeric(ctJuliaEvaluate(h, at, gradient = TRUE,
      gradient_method = "adjoint")$gradient)
    forward <- as.numeric(ctJuliaEvaluate(h, at, gradient = TRUE,
      gradient_method = "forward")$gradient)
    expect_equal(adjoint, forward, tolerance = 1e-8)
  }
})

test_that("generated counts have the dispersion's moments", {
  skip_without_julia()
  # Both generation routes, because they are separate code: the filter's and
  # the state-explicit one, the latter being what a count model's
  # `intoverstates='auto'` resolves to. A Poisson-lognormal has mean
  # `exp(mu + v/2)` and variance `mean + mean^2 (e^v - 1)`, so a generator that
  # dropped the dispersion would report a variance equal to its mean and fail
  # on the ratio alone.
  #
  # The tolerances are measured, not guessed, and they are Monte Carlo rather
  # than anything about the draw. Over twelve seeds at 15000 observations each,
  # relative error against the theoretical moments:
  #
  #   route                mean: max / typical    variance: max / typical
  #   intoverstates=TRUE      0.0107 / 0.0052        0.0358 / 0.0157
  #   intoverstates=FALSE     0.0110 / 0.0044        0.0351 / 0.0137
  #
  # So 0.05 and 0.12 are about 4.5x and 3.3x the worst seen. They were not
  # tightened when the draw became exact -- two stages, a Gaussian for the
  # dispersion and then a plain Poisson, rather than an inverted marginal --
  # because the quadrature was never what these numbers were measuring. A
  # sample variance of a heavy-tailed variable is, and 3.3x headroom on one is
  # not somewhere to economise: tightening buys no sensitivity and buys a
  # flaky test.
  sigma <- 0.5
  mu <- log(4)
  gen <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 1, manifestNames = "y", latentNames = "eta1",
    manifesttype = 3L, LAMBDA = matrix(1), DRIFT = matrix(-0.5),
    DIFFUSION = matrix(0), T0VAR = matrix(0), T0MEANS = matrix(0),
    CINT = matrix(0), MANIFESTMEANS = matrix(mu),
    MANIFESTVAR = matrix(sigma), Tpoints = 5)))
  wanted_mean <- exp(mu + sigma^2 / 2)
  wanted_var <- wanted_mean + wanted_mean^2 * (exp(sigma^2) - 1)
  for (ios in c(TRUE, FALSE)) {
    set.seed(7)
    d <- data.frame(ctGenerate(gen, n.subjects = 3000, Tpoints = 5,
      backend = "julia", intoverstates = ios))
    y <- d$y[!is.na(d$y)]
    expect_equal(mean(y), wanted_mean, tolerance = 0.05,
      info = paste("intoverstates", ios))
    expect_equal(stats::var(y), wanted_var, tolerance = 0.12,
      info = paste("intoverstates", ios))
  }
})
