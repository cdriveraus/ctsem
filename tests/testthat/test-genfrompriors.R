# ctGenerateFromPriors(), and the faults behind its intermittent failure.
#
# The \donttest example on its help page failed under `R CMD check --as-cran`
# on roughly two runs in three, with a raw stan exception -- `quad_form_sym: A
# is not symmetric. A[1,2] = nan` -- and no seed suppressed it. Two independent
# faults produced that, and each is pinned separately below, because either one
# alone leaves a hole.
#
#   1. This function fitted the model to an empty dataset so that the posterior
#      it optimised would be the prior. The priors then failed to reach that
#      fit, leaving an objective flat in every direction, and the optimiser
#      wandered until the model overflowed.
#   2. `flexlapplytext()` dropped and permuted its results at `cores > 1`, so
#      one inadmissible draw among them ended the whole call instead of being
#      skipped.
#
# The fit is gone now -- the prior is sampled directly -- so (1) cannot recur in
# that form. What replaced it is asserted here too: that the draws really are
# the prior, and that both generators turn them into data.

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

test_that("the draws are the prior, drawn rather than approximated", {
  skip_on_cran()
  # What the fit-to-empty-data route existed to produce, asserted on the thing
  # that replaced it. ctsem's raw priors are independent standard normals, so
  # this is testable as a distribution rather than as a mechanism: 200 draws,
  # checked for location, spread and independence.
  #
  # Distributional, because the failure it guards against was distributional
  # and silent. With the priors not reaching the fit, every draw came back
  # within 5e-4 of the raw origin -- the "prior predictive" was one parameter
  # vector repeated, which looks like data.
  pp <- suppressMessages(suppressWarnings(ctGenerateFromPriors(cts = ctstantestfit,
    cores = 1, nsamples = 200, parsonly = TRUE)))
  draws <- pp$stanfit$rawposterior

  expect_equal(dim(draws), c(200L, 28L))
  expect_equal(mean(draws), 0, tolerance = 0.05)
  expect_equal(stats::sd(as.numeric(draws)), 1, tolerance = 0.1)
  # Independent, not merely spread out. The largest of the 378 off-diagonal
  # correlations over 200 draws sits well inside this.
  offdiag <- abs(stats::cor(draws)[upper.tri(diag(ncol(draws)))])
  expect_lt(max(offdiag), 0.45)
})

test_that("the prior is what the model's own density says it is", {
  skip_on_cran()
  # The claim the direct draw rests on, checked against the stan program rather
  # than against a comment: with the likelihood switched off, the log density
  # over the raw parameters is the standard normal one, at the origin and away
  # from it. If a change ever gives this model a prior that is not that, this
  # fails and .ctPriorRawDraws() needs revisiting.
  #
  # `dokalman = 0` is how the likelihood is switched off; the prepared data
  # here carries the real row structure with -99 placeholders in Y, so leaving
  # it on would have the filter condition on those placeholders as if they were
  # observations. (That is what the fit-to-empty-data route arranged the long
  # way round, by handing the optimiser a dataset with nothing in it.)
  pp <- suppressMessages(suppressWarnings(ctGenerateFromPriors(cts = ctstantestfit,
    cores = 1, nsamples = 5, parsonly = TRUE)))
  npar <- ncol(pp$stanfit$rawposterior)
  prioronly <- pp$standata
  prioronly$dokalman <- 0L
  smf <- ctsem:::stan_reinitsf(pp$stanmodel, prioronly)

  expect_equal(pp$standata$priors, 1L)
  expect_equal(rstan::log_prob(smf, rep(0, npar)), npar * log(1 / sqrt(2 * pi)),
    tolerance = 1e-6)
  v <- seq(-1, 1, length.out = npar)
  expect_equal(rstan::log_prob(smf, v), sum(stats::dnorm(v, log = TRUE)),
    tolerance = 1e-6)
})

test_that("ctGenerateFromPriors() honours nsamples, and refuses `is`", {
  skip_on_cran()
  # nsamples used to be dropped on the floor: the function built an optimcontrol
  # carrying it, then overwrote the whole list two lines later, so the fit
  # always drew the stanoptimis default of 1000 however few were asked for.
  pp <- suppressMessages(suppressWarnings(ctGenerateFromPriors(cts = ctstantestfit,
    cores = 1, nsamples = 20, parsonly = TRUE)))
  expect_equal(nrow(pp$stanfit$rawposterior), 20L)
  expect_equal(nrow(pp$stanfit$transformedpars$popmeans), 20L)

  # `is` is deprecated rather than rewired, because nothing here approximates
  # anything for importance sampling to correct. It says so instead of
  # accepting quietly.
  expect_warning(
    suppressMessages(ctGenerateFromPriors(cts = ctstantestfit, cores = 1,
      nsamples = 5, parsonly = TRUE, is = TRUE)),
    regexp = 'deprecated and ignored')
})

test_that("both backends generate from the same prior draws", {
  skip_on_cran()
  # `backend` chooses the generator, not the draws. Whatever turns a parameter
  # vector into data, the shape and the completeness of the result are the
  # same, and both routes vary the dataset with the draw -- which is the point,
  # and exactly what the degenerate prior predictive did not do.
  nsamples <- 6L
  stangen <- suppressMessages(suppressWarnings(ctGenerateFromPriors(
    cts = ctstantestfit, cores = 1, nsamples = nsamples, backend = 'stan')))

  expect_equal(dim(stangen$Y)[1], nsamples)
  expect_equal(names(dimnames(stangen$Y))[1:2], c('sample', 'row'))
  expect_true(all(is.finite(stangen$Y)))
  expect_equal(length(unique(round(apply(stangen$Y, 1, mean), 9))), nsamples)
  expect_equal(dim(stangen$llrow), dim(stangen$Y)[1:2])

  skip_without_julia()
  juliagen <- suppressMessages(suppressWarnings(ctGenerateFromPriors(
    cts = ctstantestfit, cores = 1, nsamples = nsamples, backend = 'julia')))

  expect_equal(dim(juliagen$Y), dim(stangen$Y))
  expect_equal(names(dimnames(juliagen$Y))[1:2], c('sample', 'row'))
  expect_true(all(is.finite(juliagen$Y)))
  expect_equal(length(unique(round(apply(juliagen$Y, 1, mean), 9))), nsamples)
})

test_that("backend='auto' takes julia when a session is available", {
  skip_without_julia()
  # Stated as the observable consequence rather than by reading a flag: with
  # julia available, the default runs and produces what backend='julia' does.
  auto <- suppressMessages(suppressWarnings(ctGenerateFromPriors(
    cts = ctstantestfit, cores = 1, nsamples = 4)))
  expect_equal(dim(auto$Y)[c(1, 3)], c(4L, 2L))
  expect_true(all(is.finite(auto$Y)))
  expect_true(isTRUE(ctJuliaStatus()$available))
})

test_that("laplace priors are refused rather than quietly drawn as normal", {
  skip_on_cran()
  # The one thing the direct draw cannot do, and the one place it is louder
  # than what it replaced rather than merely faster. The old route drew these
  # from a gaussian approximation, so a laplaceprior parameter's "prior" draws
  # came back normal whatever laplaceprior said, and nothing reported it.
  standata <- ctstantestfit$standata
  standata$laplaceprior <- rep(1L, length(standata$laplaceprior))
  expect_error(ctsem:::.ctPriorRawDraws(standata, 28L, 5L),
    regexp = 'not a density it can draw from')
})
