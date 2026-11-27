# ctPredictVarying(): expected trajectories of hypothetical persons whose random
# effects sit at chosen multiples of their standard deviations. A random manifest
# intercept moves the observed trajectory by its sd and leaves the latent one
# alone; a random initial mean moves the latent start by its sd. Both are exact
# for intercepts and initial means, whose transforms are affine, so the checks
# compare with the fit's own random-effect sds and covariances.

.predictvarying_fit <- function() fit_cached("predictvarying_fit", {
  set.seed(31)
  data <- do.call(rbind, lapply(1:40, function(i) {
    level <- stats::rnorm(1, 0, .8)
    time <- cumsum(c(0, stats::runif(7, .4, 1.2)))
    x <- stats::rnorm(1, 0, 1)
    for (t in 2:8) {
      a <- exp(-.7 * (time[t] - time[t - 1]))
      x[t] <- a * x[t - 1] + stats::rnorm(1, 0, sqrt((1 - a^2) / 1.4))
    }
    data.frame(id = i, time = time, Y1 = x + level + stats::rnorm(8, 0, .2))
  }))
  model <- suppressWarnings(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    LAMBDA = matrix(1), MANIFESTVAR = matrix("merror"),
    MANIFESTMEANS = matrix("mm"), DRIFT = matrix("drift"),
    DIFFUSION = matrix("diff")))
  suppressMessages(ctFit(data, model, verbose = 0, cores = 1))
})

test_that("ctPredictVarying moves a person's random effects by their sds", {
  skip_without_julia()
  fit <- .predictvarying_fit()
  sds <- summary(fit)$popsd
  sds <- stats::setNames(sds$mode, rownames(sds))
  k <- suppressWarnings(ctPredictVarying(fit, vary = c("mm", "T0m_eta1"),
    plot = FALSE))
  expect_setequal(as.character(unique(k$Subject)),
    c("mm = -1", "mm = 0", "mm = 1", "T0m_eta1 = -1", "T0m_eta1 = 0", "T0m_eta1 = 1"))
  path <- function(element, subject) {
    d <- k[k$Element == element & k$Row == "Y1" & as.character(k$Subject) == subject, ]
    if (element == "etaprior") d <- k[k$Element == element & k$Row == "eta1" &
      as.character(k$Subject) == subject, ]
    d$value[order(d$Time)]
  }
  # A manifest intercept one sd up or down: the observed trajectory moves by its
  # sd at every time, the latent one not at all.
  expect_equal(path("yprior", "mm = 1") - path("yprior", "mm = 0"),
    rep(sds[["mm"]], length(path("yprior", "mm = 0"))), tolerance = 1e-6)
  expect_equal(path("yprior", "mm = 0") - path("yprior", "mm = -1"),
    rep(sds[["mm"]], length(path("yprior", "mm = 0"))), tolerance = 1e-6)
  expect_equal(path("etaprior", "mm = 1"), path("etaprior", "mm = 0"), tolerance = 1e-8)
  # An initial mean one sd up: the latent process starts higher by its sd, and
  # the difference decays.
  up <- path("etaprior", "T0m_eta1 = 1") - path("etaprior", "T0m_eta1 = 0")
  expect_equal(up[1], sds[["T0m_eta1"]], tolerance = 1e-6)
  expect_true(all(diff(up) < 0))
})

test_that("ctPredictVarying with correlated = TRUE moves the other effects by regression", {
  skip_without_julia()
  fit <- .predictvarying_fit()
  cov <- ctSummaryMatrices(fit)$estimate$popcov
  k <- suppressWarnings(ctPredictVarying(fit, vary = "mm", correlated = TRUE,
    plot = FALSE))
  latent <- function(subject) {
    d <- k[k$Element == "etaprior" & as.character(k$Subject) == subject, ]
    d$value[order(d$Time)][1]
  }
  # The initial mean takes its expected value given a manifest intercept one sd
  # up: its covariance with the intercept over the intercept's sd.
  expect_equal(latent("mm = 1") - latent("mm = 0"),
    cov["T0m_eta1", "mm"] / sqrt(cov["mm", "mm"]), tolerance = 1e-6)
})

test_that("ctPredictVarying draws a plot for each name and refuses what it cannot vary", {
  skip_without_julia()
  fit <- .predictvarying_fit()
  plots <- suppressWarnings(ctPredictVarying(fit, sdValues = c(-2, 0, 2)))
  expect_setequal(names(plots$Process$Observed), c("T0m_eta1", "mm"))
  expect_true(inherits(plots$Process$Latent$mm, "ggplot"))
  expect_error(ctPredictVarying(fit, vary = "drift", plot = FALSE),
    "neither a time independent predictor nor a parameter")
  expect_error(ctPredictTIP(fit, plot = FALSE), "no time independent predictors")
})
