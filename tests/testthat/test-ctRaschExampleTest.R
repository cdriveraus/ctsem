# A Rasch-style measurement model: seven binary indicators loading on one
# latent process, item difficulties as free manifest means with the first fixed
# at zero for identification, and a per-subject continuous intercept.
#
# All of it runs on the julia backend. Nothing here is kept on stan, and the
# reasons are specific rather than a preference:
#
#   - The measurement update is the point of the file. Stan moment-matches the
#     Bernoulli observation to a Gaussian; julia integrates it. Comparing the
#     two on a binary model measures the linearisation, not either estimator,
#     which is why the stan-to-stan comparison this file used to make needed a
#     tolerance of 0.2 to pass -- about one and a half posterior standard
#     deviations, loose enough that it could not have failed for a real reason.
#   - The stan binary path is still exercised end to end, in
#     test-binary-path.R, which is where a regression in it would show.
#   - A julia fit needs no stan compile, which is most of what this file used
#     to cost.
#
# What replaced the old comparison: augmented ML against Laplace ML, two
# estimators of the same objective that differ in how the random effects are
# handled, compared on the raw parameter vector so no Monte Carlo enters. They
# agree to 0.007 -- far tighter than the check they replace.

# The generating process, in base R so that the data move only with R's own
# random stream: ctGenerate's draw order changes with commits that have nothing
# to do with this file, and its data here once did. One latent process with
# drift -0.3 and diffusion sd 0.3, integrated exactly over intervals
# exp(N(0, 0.2^2)) long, 20 occasions a subject; each subject's intercept from
# N(0.1, 0.3^2) and its initial state from N(-0.3, 0.5^2). Seven binary
# indicators through the logistic link the fit uses, P(y = 1) =
# plogis(eta + difficulty), the first difficulty zero. The "not parked at zero"
# assertion below is what would show data generated without the individual
# variation, which this file once spent a long time fitting.
.rasch_data <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    set.seed(1234)
    nsubjects <- 20
    nt <- 20
    difficulty <- c(m2 = .5, m3 = .5, m4 = .5, m5 = -.5, m6 = -.5, m7 = -.5)
    drift <- -.3
    diffusion <- .3
    intercept <- stats::rnorm(nsubjects, .1, .3)
    d <- do.call(rbind, lapply(seq_len(nsubjects), function(i) {
      dt <- exp(stats::rnorm(nt - 1, 0, .2))
      eta <- numeric(nt)
      eta[1] <- stats::rnorm(1, -.3, .5)
      for (t in 2:nt) {
        e <- exp(drift * dt[t - 1])
        eta[t] <- e * eta[t - 1] + (e - 1) / drift * intercept[i] +
          stats::rnorm(1, 0, diffusion * sqrt((e^2 - 1) / (2 * drift)))
      }
      p <- stats::plogis(outer(eta, c(0, difficulty), "+"))
      y <- matrix(as.numeric(stats::runif(length(p)) < p), nt)
      colnames(y) <- paste0("Y", seq_len(ncol(y)))
      data.frame(id = i, time = c(0, cumsum(dt)), y)
    }))
    cache <<- list(
      data = d,
      # The realised spread, not the nominal 0.3: 20 draws of an sd 0.3 normal
      # have a sample sd of their own, and it is that the fit can recover.
      intercept_sd = stats::sd(intercept),
      # Item difficulties, in the order the fitted model names them m2..m7.
      difficulty = difficulty,
      drift = drift, diffusion = diffusion)
    cache
  }
})

# The first manifest mean is fixed at zero; the rest are the item difficulties,
# with an identity transform, so a raw parameter m2..m7 *is* a difficulty and
# the fits below can be compared to each other without going through draws.
.rasch_model <- function(type = "ct") {
  n.manifest <- 7
  m <- suppressMessages(ctModel(n.latent = 1, n.manifest = n.manifest,
    MANIFESTMEANS = c(0, paste0('m', 2:n.manifest, '|param|FALSE')),
    LAMBDA = rep(1, n.manifest),
    DIFFUSION = 'diff|log1p_exp(2*param)',
    T0MEANS = 't0m|param|TRUE|.2',
    CINT = 'b|param|TRUE|1',
    type = type))
  m$manifesttype[] <- 1L   # binary
  m
}

# Three fits, memoised: two of the four tests below compare fits to each other,
# and the julia engine's first call in a session pays its compilation. Refitting
# per test would roughly treble the file.
#
# `priors = TRUE` throughout. The model's own labels ask for priors -- 'param'
# with no distribution named is N(0,1) -- and they were inert here for years
# because every optimised fit passed `priors = FALSE`. They are not decoration:
# at `priors = FALSE` the discrete-time fit is free to trade its autoregression
# against its innovation variance along a nearly flat ridge, and from a zero
# start it walks to an autoregression of 0.97 with a diffusion of 7.5e5, for
# 0.8 of a nat.
.rasch_fit <- local({
  cache <- list()
  function(what) {
    if (!is.null(cache[[what]])) return(cache[[what]])
    d <- .rasch_data()$data
    args <- switch(what,
      # The fit everything else is read against, and the only one that needs
      # draws: its intervals are what the recovery checks are stated in.
      ct = list(model = .rasch_model("ct"), intoverpop = "augmented"),
      # Laplace integrates the random effects instead of maximising over them
      # jointly. Compared on the raw vector, so 20 draws is all it needs -- and
      # the draws are most of a Laplace fit's cost, because each one redoes the
      # inner solves.
      laplace = list(model = .rasch_model("ct"), intoverpop = "laplace",
        optimcontrol = list(finishsamples = 20)),
      # Same random-effect treatment as `ct`, so the pair differs in exactly
      # one thing: the discretization.
      dt = list(model = .rasch_model("dt"), intoverpop = "augmented",
        optimcontrol = list(finishsamples = 20)))
    fit <- suppressWarnings(suppressMessages(do.call(ctFit, c(
      list(datalong = d, backend = "julia", cores = 1, verbose = 0,
        intoverstates = TRUE, priors = TRUE, optimize = TRUE), args))))
    cache[[what]] <<- fit
    fit
  }
})

.rasch_raw <- function(fit) {
  raw <- fit$estimate$raw
  stats::setNames(as.numeric(raw),
    ctsem:::.ctBackendRawParameterNames(fit, length(raw)))
}

test_that("a Rasch model recovers the item difficulties that generated it", {
  skip_on_cran()
  skip_without_julia()
  truth <- .rasch_data()$difficulty
  est <- suppressWarnings(suppressMessages(summary(.rasch_fit("ct"))))$popmeans
  items <- names(truth)

  expect_true(all(items %in% rownames(est)))
  # Intervals rather than point estimates. Seven binary indicators at 20
  # occasions leave a difficulty with a posterior sd around 0.15, so a point
  # estimate is allowed to be off by that much; what must hold is that the
  # model's own uncertainty covers the truth. The measured margins are
  # comfortable, the tightest being m3 at 0.276 [-0.009, 0.555] against a
  # truth of 0.5.
  for (it in items) {
    expect_lt(est[it, "2.5%"], truth[[it]], label = paste0(it, " lower bound"))
    expect_gt(est[it, "97.5%"], truth[[it]], label = paste0(it, " upper bound"))
  }
  # And the easy/hard split is recovered with no tolerance at all: the three
  # items generated at +0.5 come back positive and the three at -0.5 negative.
  # An estimator that had lost the item structure -- rather than merely being
  # imprecise about it -- would fail this while still covering above.
  expect_true(all(est[c("m2", "m3", "m4"), "mean"] > 0))
  expect_true(all(est[c("m5", "m6", "m7"), "mean"] < 0))
})

test_that("the between-subject spread in the intercept is recovered, not parked at zero", {
  skip_on_cran()
  skip_without_julia()
  truth <- .rasch_data()
  s <- suppressWarnings(suppressMessages(summary(.rasch_fit("ct"))))

  expect_true("b" %in% rownames(s$popsd))
  # The realised sd of the 20 generated intercepts is 0.304; the fit reports
  # 0.35 [0.16, 0.67].
  expect_lt(s$popsd["b", "2.5%"], truth$intercept_sd)
  expect_gt(s$popsd["b", "97.5%"], truth$intercept_sd)
  # The assertion that would have caught the generating bug, which coverage
  # alone would not: it is the distance from zero, rather than the interval's
  # width, that says the random effect is there at all. When every subject is
  # generated from the same intercept this fit returns 0.071 [0.003, 0.271],
  # so 0.05 sits between the two measured lower bounds, 0.003 and 0.161.
  expect_gt(s$popsd["b", "2.5%"], 0.05)

  # The dynamics, on the same terms. Wide, because seven binary indicators
  # carry much less than one continuous one: drift -0.328 [-0.617, -0.138]
  # against -0.3, diffusion 0.230 [0.026, 0.448] against 0.3.
  expect_lt(s$popmeans["drift_eta1", "2.5%"], truth$drift)
  expect_gt(s$popmeans["drift_eta1", "97.5%"], truth$drift)
  expect_lt(s$popmeans["diff", "2.5%"], truth$diffusion)
  expect_gt(s$popmeans["diff", "97.5%"], truth$diffusion)
})

test_that("the augmented and Laplace estimators find the same optimum", {
  skip_on_cran()
  skip_without_julia()
  # The heir of this file's one real check -- two estimators, same data -- moved
  # onto the backend that integrates the measurement rather than linearising it,
  # and onto the raw vector, where no draw noise stands between them.
  #
  # They are genuinely different objectives: `intoverpop = 'augmented'`
  # integrates the random effects through the filter's linearised measurement
  # update, `'laplace'` by a Laplace approximation at each subject's mode. Agreement is
  # therefore a result and not an identity, and it is not exact: the largest
  # gap is 0.007, on the drift. The item difficulties agree to 0.0003.
  augmented <- .rasch_raw(.rasch_fit("ct"))
  laplace <- .rasch_raw(.rasch_fit("laplace"))

  expect_identical(names(augmented), names(laplace))
  expect_lt(max(abs(augmented - laplace)), 0.1)
  items <- paste0("m", 2:7)
  expect_lt(max(abs(augmented[items] - laplace[items])), 0.01)
  # Two optima of two different objectives, so the likelihoods are close rather
  # than equal: -1720.70 augmented, -1720.58 as the Laplace fit reports it.
  expect_lt(abs(.rasch_fit("ct")$estimate$loglik -
      .rasch_fit("laplace")$estimate$loglik), 2)
})

test_that("the discrete-time twin sees the same measurement model", {
  skip_on_cran()
  skip_without_julia()
  # What `rod`/`sod` were presumably for. They were fitted and summarised here
  # for years and never asserted on, so a third of the file's runtime produced
  # nothing.
  #
  # The correspondence that actually holds is about the measurement model, not
  # the dynamics: an item difficulty is identified by the Bernoulli link and
  # the loadings, neither of which knows how time was discretised, so the two
  # fits must agree on it. They agree to 0.0005.
  continuous <- .rasch_raw(.rasch_fit("ct"))
  discrete <- .rasch_raw(.rasch_fit("dt"))
  items <- paste0("m", 2:7)

  # Same parameters in the same order -- the discretization changes what DRIFT
  # means, not how many free cells the model has.
  expect_identical(names(continuous), names(discrete))
  expect_lt(max(abs(continuous[items] - discrete[items])), 0.05)

  # The dynamics parameter is *not* shared, and asserting that is what keeps
  # the check above from being vacuous. In discrete time DRIFT is already the
  # one-step transition, so it must come back as an autoregression in (0, 1):
  # about 0.75 here, against a continuous rate of -0.328.
  s_ct <- suppressWarnings(suppressMessages(summary(.rasch_fit("ct"))))
  s_dt <- suppressWarnings(suppressMessages(summary(.rasch_fit("dt"))))
  expect_lt(s_ct$popmeans["drift_eta1", "mean"], 0)
  expect_gt(s_dt$popmeans["drift_eta1", "mean"], 0)
  expect_lt(s_dt$popmeans["drift_eta1", "mean"], 1)

  # Deliberately not asserted: that the discrete transition equals
  # exp(drift * dt). It does not reliably. The observed intervals average 1.01
  # but vary (sd 0.2 on the log scale), and the discrete likelihood is nearly flat along a
  # ridge trading the autoregression against the innovation variance -- over the
  # configurations measured here the discrete fit landed anywhere from 0.69 to
  # 0.97 while the continuous prediction stayed near 0.76, for differences in
  # log likelihood under one nat. Pinning that number would be pinning which end
  # of a ridge the optimizer stopped at.
  #
  # What is stable is that the two discretizations describe this data about
  # equally well: -1720.70 continuous against -1720.98 discrete.
  expect_lt(abs(.rasch_fit("ct")$estimate$loglik -
      .rasch_fit("dt")$estimate$loglik), 5)
})
