# The population summaries -- `popmeans`, `popsd`, `popcov` -- are the only
# place the stan program reports a random effect on the scale a user reads.
# Everything else about a fit can be right while they are wrong, because they
# live in generated quantities and never touch the likelihood: the fit, the log
# likelihood and every raw-scale quantity (`rawpopmeans`, `rawpopsd`,
# `rawpopcorr`) are unaffected by anything asserted here. That is exactly what
# makes the failure quiet, so it gets its own fast test rather than being left
# to the multi-minute simulation files that first caught it.
#
# The specific thing at stake: with `intoverpop`, an indvarying parameter is
# carried by an extra latent state whose T0MEANS cell holds the *raw* value,
# and the parameter's real transform lives in the cell it feeds -- `10*state[3]`
# in CINT, or `-(1e-6 + 2*log1p_exp(-2*state[2]))` in DRIFT. The generated
# quantities have to follow that reference (matsetup's `stateref` column) to
# find the transform. When they do not, every such parameter is summarised
# untransformed: a CINT population sd ten times too small, a DRIFT one on the
# unconstrained scale entirely, while correlations -- being scale free -- still
# look perfect.
#
# No fit is needed, and none is done. `ctFit(fit=FALSE)` plus
# `rstan::constrain_pars` evaluates the generated quantities at raw values we
# choose, so the expected numbers are arithmetic rather than estimates.

skip_on_cran()
{  # body of the guard this replaced; indentation unchanged

  library(ctsem)
  library(testthat)

  context("popsummaryscale")

  test_that("population summaries use the transform of the cell a random effect feeds", {
    skip_on_cran()
    skip_if_not_installed("rstan")

    # cint is indvarying by default and reaches CINT through `10*state[3]`;
    # drift is made indvarying and reaches DRIFT through the default
    # `-(1e-6 + 2*log1p_exp(-2*state[2]))`. One linear cell, one not.
    model <- ctModel(type = 'ct', silent = TRUE, LAMBDA = matrix(1),
      DRIFT = matrix(ctParSpec('drift', indvarying = TRUE)),
      DIFFUSION = matrix(.5), MANIFESTVAR = matrix(.5),
      MANIFESTMEANS = matrix(0), CINT = matrix('cint'),
      T0MEANS = matrix(0), T0VAR = matrix(1))
    datalong <- data.frame(id = rep(1:3, each = 3), time = rep(c(0, .5, 1.5), 3),
      Y1 = c(0, .1, .2, .1, 0, -.1, .2, .3, .1))

    spec <- suppressMessages(ctFit(datalong, model, backend = 'stan',
      fit = FALSE, priors = FALSE))
    # Asserted rather than branched on: this shape must reuse the precompiled
    # binary. If it ever stops doing so the test would silently start
    # compiling a stan program, which is minutes, not seconds.
    expect_equal(as.integer(spec$standata$recompile), 0L)
    expect_equal(as.integer(spec$standata$nindvarying), 2L)

    sf <- ctsem:::stan_reinitsf(ctsem:::stanmodels$ctsm, spec$standata)

    # raw vector layout (see .ctBackendLaplacePriorSpec in ctBackendUncertainty.R):
    # nparams rawpopmeans, then nindvarying rawpopsdbase, then the correlation
    # coordinates. rawpopsd = log1p_exp(2*rawpopsdbase-1), inverted here so the
    # population sds are round numbers; the single correlation coordinate is set
    # to 0, which is a correlation of zero, so the two effects stay independent
    # and each expected value is a one dimensional calculation.
    rawmeans <- c(drift = 1.2, cint = .3)
    rawsd <- c(drift = .5, cint = .3)
    upars <- c(rawmeans, (log(exp(rawsd) - 1) + 1) / 2, 0)

    cp <- rstan::constrain_pars(sf, upars)
    expect_equal(as.numeric(cp$rawpopsd), as.numeric(rawsd), tolerance = 1e-8)

    # Means are exact: no simulation is involved, so a population mean must
    # equal the matrix cell it is the population mean of.
    # 1e-4 rather than exact: the two are computed by different routes (this one
    # by tform on the raw value, pop_DRIFT by the model block from the carrier
    # state) and differ in the last few digits. Untransformed, popmeans[1] would
    # be 1.2.
    expect_equal(as.numeric(cp$popmeans[1]), as.numeric(cp$pop_DRIFT[1, 1]), tolerance = 1e-4)
    expect_equal(as.numeric(cp$popmeans[2]), as.numeric(cp$pop_CINT[1, 1]), tolerance = 1e-8)
    expect_equal(as.numeric(cp$popmeans[2]), 10 * rawmeans[['cint']], tolerance = 1e-8)

    # Standard deviations are a 1000 draw monte carlo inside the stan program
    # (`popcovn`), redrawn on every call, so average a few and keep the
    # tolerances loose. They are still nowhere near loose enough to admit the
    # untransformed answers: 0.3 instead of 3 for cint, and 0.59 instead of
    # 0.27 for drift.
    reps <- t(replicate(5, as.numeric(rstan::constrain_pars(sf, upars)$popsd)))
    popsd <- colMeans(reps)

    # Linear cell: the population sd is the raw sd times the cell's multiplier.
    expect_equal(popsd[2], 10 * rawsd[['cint']], tolerance = .15)

    # Nonlinear cell: the population sd is the spread of the transform over the
    # population, and the population is centred at rawpopmeans. Centring it
    # anywhere else -- at rawpopcovchol * rawpopmeans, say -- reads the
    # transform's slope at the wrong point and gives 0.59 here.
    drifttf <- function(param) -(1e-6 + 2 * log1p(exp(-(2 * param))))
    set.seed(1)
    z <- rnorm(2e5)
    expect_equal(popsd[1], sd(drifttf(rawmeans[['drift']] + rawsd[['drift']] * z)),
      tolerance = .2)
  })
}
