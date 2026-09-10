skip_on_cran()
skip_on_32bit()
{  # body of the guard this replaced; indentation unchanged
  library(ctsem)
  library(testthat)

  context("ctBinaryBinaryMix")

  # Two latents observed only through binary indicators, asymmetrically: eta1
  # through one, eta2 through ten. That asymmetry is the whole design -- it is
  # the heavily indicated latent whose process noise a linearised measurement
  # update loses, so a fit that recovers eta1 and not eta2 is diagnostic rather
  # than merely noisy.
  #
  # THE GENERATING MODEL IS WRITTEN OUT IN FULL, AND MUST STAY THAT WAY.
  # It used not to be: MANIFESTVAR, T0VAR, T0MEANS and MANIFESTMEANS were left
  # free and `ctGenerate()` filled them from `.ctGenerateDefaults()`. Those
  # defaults changed in 4a27c253 (MANIFESTVAR diagonal 0.5 -> 0, T0VAR 1 ->
  # 1e-6) and the data this file fits changed with them, silently, under
  # assertions written against the old data and reading as though the seed
  # fixed everything. Measured: the same code at the two sets of defaults gives
  # DIFFUSIONcov[2,2] of 0.098 [0.081, 0.118] and 0.148 [0.123, 0.180], against
  # a true 0.16. Nothing about the estimator moved between those two numbers.
  #
  # The old defaults were also the wrong data to be testing on. A generating
  # MANIFESTVAR of 0.5 puts residual variance on a manifest that the fitted
  # model cannot have -- ctFit() fixes a binary indicator's MANIFESTVAR to zero
  # by construction -- so the low process-noise estimate was measuring that
  # misspecification as much as the measurement update. Generating with
  # MANIFESTVAR zero leaves the measurement update as the only approximation
  # under study, which is what this file is for.

  .bbmix_gen <- function(seed, n.subjects = 50, Tpoints = 200, nind = 10){
    invlog <- function(x) exp(x)/(1 + exp(x))
    set.seed(seed)
    gm <- ctModel(DRIFT = c(-.2, .2,
      0, -.1),
      DIFFUSION = c(.3, 0,
        0, .4),
      CINT = c(.1, .1),
      LAMBDA = diag(1, 2),
      MANIFESTVAR = diag(0, 2),
      MANIFESTMEANS = matrix(0, 2, 1),
      T0VAR = diag(1, 2),
      T0MEANS = matrix(0, 2, 1),
      n.latent = 2, n.manifest = 2, Tpoints = Tpoints)
    d <- ctGenerate(gm, n.subjects = n.subjects, logdtsd = .1, dtmean = .1,
      burnin = 20)
    d[, 'Y1'] <- rbinom(nrow(d), size = 1, prob = invlog(d[, 'Y1']))
    d <- data.frame(d)
    for(i in seq_len(nind)){
      d[[paste0('b', i)]] <- rbinom(nrow(d), size = 1, prob = invlog(d[, 'Y2']))
    }
    # The generating matrices, read the way the current model object holds
    # them: `pars`, not top-level fields. `gm$DIFFUSION` was NULL for years
    # here and nobody saw it, because the file was named `ctBinaryGaussianMix.R`
    # and testthat only runs `test-*`.
    gmn <- ctModelMatrices(ctsem:::ctModeltoNumeric(gm))
    gmn$DIFFUSIONcov <- tcrossprod(gmn$DIFFUSION)
    list(d = d, gmn = gmn, nind = nind)
  }

  .bbmix_model <- function(nind = 10){
    MANIFESTVAR <- diag(c(1, rep(0, nind)), nind + 1)
    MANIFESTVAR[1] <- 'mvar1'
    m <- ctModel(type = 'ct',
      manifestNames = c('Y1', paste0('b', seq_len(nind))),
      LAMBDA = rbind(diag(1, 2), cbind(rep(0, nind - 1), rep(1, nind - 1))),
      MANIFESTMEANS = 0,
      MANIFESTVAR = MANIFESTVAR,
      CINT = c('CINT1', 'cint2'))
    m$manifesttype[seq_len(nind + 1)] <- 1 #set type to binary
    m$pars$indvarying <- FALSE
    m
  }

  # Recovery, stated as coverage of the generating value by the 95% interval,
  # for every population quantity this design identifies. Coverage rather than a
  # tolerance because the interval is the fit's own statement of what it knows;
  # a tolerance chosen to fit the answer would pass whatever came out.
  .bbmix_expect_recovery <- function(f, gmn, label){
    low <- ctSummaryMatrices(f, calcfuncargs = list(probs = .025))
    up  <- ctSummaryMatrices(f, calcfuncargs = list(probs = .975))
    cells <- list(c('DRIFT', 1, 1), c('DRIFT', 2, 1), c('DRIFT', 1, 2),
      c('DRIFT', 2, 2), c('CINT', 1, 1), c('CINT', 2, 1),
      c('DIFFUSIONcov', 1, 1), c('DIFFUSIONcov', 2, 2))
    for(cell in cells){
      nm <- cell[1]; i <- as.integer(cell[2]); j <- as.integer(cell[3])
      expect_true(low[[nm]][i, j] < gmn[[nm]][i, j] &&
          up[[nm]][i, j] > gmn[[nm]][i, j],
        info = sprintf('%s: %s[%d,%d] true %.4f, interval [%.4f, %.4f]',
          label, nm, i, j, gmn[[nm]][i, j], low[[nm]][i, j], up[[nm]][i, j]))
    }
    invisible(NULL)
  }

  test_that("ctBinaryBinaryMix1", {
    # Recovery, not a characterisation of bias. This file used to assert the
    # bias instead -- that DIFFUSIONcov[2,2] was low enough to miss its own
    # interval, and that a truly-zero cross-effect came back significant -- and
    # that was right for the data it was then generating. It is not right for
    # this data: with the generating model written out above, the stan fit
    # covers all eight quantities, checked at seeds 1234, 1 and 2 (24 of 24).
    #
    # This is NOT evidence that the linearised update improved; nothing in the
    # measurement path of `inst/stan/ctsm.stan` has changed. It is a design in
    # which the approximation does not bite: 10000 observations, densely
    # sampled, correctly specified. What the approximation costs is recorded
    # where it was measured -- the RMSE table in `R/ctFit.R` above the binary
    # message, and `test-julia-binary.R`. Nor is that cost something a
    # one-dataset assertion can pin: on the old data the spurious cross-effect
    # this file used to assert appeared at seed 1234 and at neither of seeds 1
    # and 2, so that expectation was measuring a seed.
    #
    # The stan half stays here rather than moving to `dev/STAN-DEPRECATION.md`:
    # stan is still the default backend, and this is the only test that fits
    # eleven binary indicators through it. It runs under CTSEM_TEST_STAN now
    # rather than on every run: it is about half the file's 295 s, on the
    # backend being deprecated, for the same eight recovery claims the julia
    # block below makes. The two already share `.bbmix_expect_recovery()`, so
    # this is the same test on both backends whenever it runs; the flag
    # decides only how often the expensive half is paid for.
    skip_if_not(.ctsem_test_stan(), "CTSEM_TEST_STAN is not set.")
    gen <- .bbmix_gen(1234)
    f <- ctFit(datalong = gen$d, model = .bbmix_model(gen$nind), cores = 2,
      plot = 10)
    .bbmix_expect_recovery(f, gen$gmn, 'stan')
  })

  test_that("ctBinaryBinaryMixJulia", {
    skip_without_julia()
    # The same design on the backend that integrates the Bernoulli observation
    # against the predicted state instead of moment-matching it to a Gaussian.
    # This shape -- two latents, one of them behind ten binary indicators -- is
    # the one the linearisation was measured to lose, and it was not covered on
    # julia anywhere: `test-julia-binary.R` fits a single latent, and its mixed
    # test pairs binary indicators with a continuous one. Checked at seeds
    # 1234, 1 and 2, as the stan half was: 24 of 24 covered on each backend.
    #
    # Full size deliberately. Everything smaller was tried and eta1, which has
    # a single binary indicator, stops being identified: at 25 subjects x 200
    # timepoints its DIFFUSIONcov missed at one seed of two, and at 25 x 60 the
    # interval came back [0, 5364]. About four minutes locally, which is what
    # the design costs; `optimcontrol = list(estonly = TRUE)` saves none of it,
    # since the optimisation and not the uncertainty phase is the cost.
    gen <- .bbmix_gen(1234)
    f <- ctFit(datalong = gen$d, model = .bbmix_model(gen$nind), cores = 1,
      backend = 'julia')
    .bbmix_expect_recovery(f, gen$gmn, 'julia')
  })

}
