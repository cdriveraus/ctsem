# ctGenerateFromPriors(), and the two defects behind its intermittent failure.
#
# The \donttest example on its help page failed under `R CMD check --as-cran`
# on roughly two runs in three, with a raw stan exception -- `quad_form_sym: A
# is not symmetric. A[1,2] = nan` -- and no seed suppressed it. Two independent
# faults produced that, and each is pinned separately below, because either one
# alone leaves a hole.
#
#   1. The priors never reached the fit, so the objective was flat: no data and
#      no priors. The optimizer then had no gradient to hold it anywhere, and
#      the draws it eventually produced could be anywhere.
#   2. `flexlapplytext()` dropped and permuted its results at `cores > 1`, so
#      one inadmissible draw among them ended the whole call instead of being
#      skipped.
#
# The first is the silent one and matters on its own: a flat objective gives a
# "prior predictive" that is one parameter vector repeated, which looks like
# data and is not.

test_that("flexlapplytext() returns results in input order and keeps NULLs", {
  skip_on_cran()
  # The mechanism, with no stan in the way. `stan_constrainsamples()` reads
  # both properties of this list: it locates the first admissible sample by
  # position, and it uses NULL to mean "this draw was rejected". At cores > 1
  # the results used to come back concatenated in worker order with every NULL
  # silently removed by unlist(), so a rejected first draw was invisible and
  # the position no longer named the sample it came from.
  cl <- ctsem:::makeClusterID(2)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  ctsem:::clusterIDeval(cl, list('tfun <- function(i) if(i %% 3 == 0) NULL else i * 10'))

  n <- 10
  out <- ctsem:::flexlapplytext(cl, 1:n, 'tfun', cores = 2)
  expected <- lapply(1:n, function(i) if(i %% 3 == 0) NULL else i * 10)

  expect_equal(length(out), n)
  expect_equal(out, expected)
  # Stated separately from the equality above: this is the property the caller
  # actually depends on, and an off-by-one in the reassembly would satisfy
  # neither but is easiest to read here.
  expect_equal(which(vapply(out, is.null, logical(1))), c(3L, 6L, 9L))

  # cores = 1 was always correct; it is the reference the parallel branch has
  # to match, so say so rather than assuming it. That branch evaluates the
  # function name in the caller's frame rather than in a worker, so it needs
  # the definition here.
  tfun <- function(i) if(i %% 3 == 0) NULL else i * 10
  expect_equal(ctsem:::flexlapplytext(cl, 1:n, 'tfun', cores = 1), expected)
})

test_that("stan_constrainsamples() skips an inadmissible first draw at cores>1", {
  skip_on_cran()
  # The user-facing consequence. A draw the stan program refuses -- here a raw
  # vector far enough out that the filter covariance stops being positive
  # definite -- must be dropped and reported, not thrown. Row 1 in particular,
  # because that is the row the skeleton was built from unprotected: with the
  # NULLs dropped, `which(!nulls)[1]` was always 1 no matter which draw failed.
  fit <- ctstantestfit
  npar <- length(fit$stanfit$rawest)
  standata <- fit$standata
  standata$savescores <- 0L
  standata$savesubjectmatrices <- 0L

  smf <- ctsem:::stan_reinitsf(fit$stanmodel, standata)
  bad <- rep(20, npar)
  skip_if_not(inherits(try(rstan::constrain_pars(smf, upars = bad), silent = TRUE),
    'try-error'), 'no inadmissible raw vector available for this model')

  samples <- rbind(bad, matrix(fit$stanfit$rawest, nrow = 4, ncol = npar, byrow = TRUE))
  out <- suppressMessages(suppressWarnings(ctsem:::stan_constrainsamples(
    sm = fit$stanmodel, standata = standata, samples = samples, cores = 2,
    savescores = FALSE, savesubjectmatrices = FALSE, dokalman = TRUE,
    onlyfirstrow = FALSE, pcovn = FALSE)))

  # Four of the five survive, and the returned draws are the admissible ones.
  expect_equal(dim(out$rawpopmeans)[1], 4L)
})

test_that("ctGenerateFromPriors() fits the empty dataset with priors on", {
  skip_on_cran()
  # `priors` was read only from `cts$args$resolved$priors`, a field that fits
  # made before it existed -- ctstantestfit among them -- do not carry. NULL
  # there meant `args$priors <- NULL` *removed* the element rather than setting
  # it, so ctFit() fell back to its own default of FALSE and the guard for
  # exactly this case never fired.
  #
  # Measured on the unfixed code: log density identically 0 at every raw
  # vector tried, a Hessian with no direction of positive curvature, a repaired
  # covariance of 1e-8 * I, and all 1000 draws within 5e-4 of the raw origin.
  # Every generated dataset came from the same parameter vector.
  expect_null(ctstantestfit$args$resolved)

  pp <- suppressMessages(suppressWarnings(ctGenerateFromPriors(cts = ctstantestfit,
    cores = 1, nsamples = 20, parsonly = TRUE)))

  expect_equal(pp$standata$priors, 1L)
  # The raw priors are normal(0,1), so the draws carry that spread. The point
  # is the contrast with 5e-4, not the exact number.
  expect_gt(stats::sd(as.numeric(pp$stanfit$rawposterior)), 0.5)
  # And the objective is the prior, not a flat surface.
  smf <- ctsem:::stan_reinitsf(pp$stanmodel, pp$standata)
  npar <- length(pp$stanfit$rawest)
  expect_lt(rstan::log_prob(smf, rep(1, npar)), rstan::log_prob(smf, rep(0, npar)))
})

test_that("ctGenerateFromPriors() honours nsamples and is", {
  skip_on_cran()
  # Both used to be dropped on the floor. The function built an optimcontrol
  # carrying `is` and finishsamples, then overwrote the whole list two lines
  # later, so the fit always drew the stanoptimis default of 1000 no matter what
  # nsamples said, and `is` did nothing. Wiring the old `optimcontrol$is`
  # through would not have worked either -- ctFit() refuses that name outright,
  # so the call would have stopped rather than importance sampled.
  hess <- suppressMessages(suppressWarnings(ctGenerateFromPriors(cts = ctstantestfit,
    cores = 1, nsamples = 20, parsonly = TRUE)))
  expect_equal(nrow(hess$stanfit$rawposterior), 20L)
  expect_equal(hess$stanfit$uncertainty$settings$method, 'hessian')

  isfit <- suppressMessages(suppressWarnings(ctGenerateFromPriors(cts = ctstantestfit,
    cores = 1, nsamples = 20, parsonly = TRUE, is = TRUE)))
  expect_equal(isfit$stanfit$uncertainty$settings$method, 'is')
  expect_equal(nrow(isfit$stanfit$rawposterior), 20L)

  # And it is the same prior either way, which is the point of the note on the
  # argument: with every raw prior normal(0,1) and no data, the target is
  # exactly gaussian, so there is nothing for the importance weights to
  # correct. Measured sd 1.006 against 1.009.
  expect_equal(stats::sd(as.numeric(isfit$stanfit$rawposterior)), 1, tolerance = .15)
  expect_equal(stats::sd(as.numeric(hess$stanfit$rawposterior)), 1, tolerance = .15)
})
