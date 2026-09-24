# ctVarianceDecomposition -----------------------------------------------------
#
# The checks that matter here are against independently derived answers, not
# against the other backend: a decomposition built out of the same matrices the
# filter uses would agree with itself while both were wrong. So the stationary
# case is checked against the closed form written out in the test
# (asymDIFFUSIONcov = -DIFFUSIONcov / (2 DRIFT) for one process), the
# non-stationary case against a transient computed from expm() here, and the
# between person variance against the random effect spread the model was
# generated with.

test_that('Gauss-Hermite nodes integrate a standard normal', {
  gh <- ctsem:::.ctVarDecompGaussHermite(21L)
  expect_equal(sum(gh$weight), 1)
  expect_equal(sum(gh$weight * gh$node), 0)
  expect_equal(sum(gh$weight * gh$node^2), 1)
  expect_equal(sum(gh$weight * gh$node^4), 3)
  expect_error(ctsem:::.ctVarDecompGaussHermite(2L), 'at least 3')
})

test_that('one interval step matches the closed forms it generalises', {
  drift <- matrix(c(-0.5, 0, 0.2, -0.3), 2, 2)
  diffusioncov <- matrix(c(1, 0.2, 0.2, 0.6), 2, 2)
  cint <- matrix(c(0.7, -0.4), 2, 1)
  dt <- 1.3
  step <- ctsem:::.ctVarDecompStep(drift, diffusioncov, cint, dt, TRUE)
  expect_equal(step$transition, as.matrix(expm::expm(drift * dt)))
  # The block exponential intercept is solve(DRIFT) %*% (expm(DRIFT dt) - I) %*%
  # CINT wherever DRIFT is invertible, which is ctsem's own discreteCINT.
  expect_equal(step$intercept,
    solve(drift) %*% (as.matrix(expm::expm(drift * dt)) - diag(2)) %*% cint)
  # And the innovation over an infinite interval is the asymptotic covariance.
  asym <- matrix(solve(kronecker(diag(2), drift) + kronecker(drift, diag(2)),
    -as.vector(diffusioncov)), 2, 2)
  expect_equal(ctsem:::.ctVarDecompStep(drift, diffusioncov, cint, 200, TRUE)$innovation,
    asym, tolerance = 1e-6)
})

test_that('a discrete time step accumulates the intercept over its steps', {
  drift <- matrix(c(0.5, 0, 0.1, 0.4), 2, 2)
  cint <- matrix(c(1, -2), 2, 1)
  step <- ctsem:::.ctVarDecompStep(drift, diag(2) * 0.5, cint, 3, FALSE)
  expect_equal(step$transition, drift %*% drift %*% drift)
  expect_equal(step$intercept, cint + drift %*% cint + drift %*% drift %*% cint)
})


# Fits ------------------------------------------------------------------------

skip_without_julia()

# One process, stationary start, a random intercept. Everything this file needs
# a closed form for is a scalar here.
set.seed(11)
truedrift <- -0.5
truediffusion <- 1
truemanifestvar <- 0.4
traitsd <- 0.3
nsub <- 25
tpoints <- 20

generating <- ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
  LAMBDA = matrix(1), DRIFT = matrix(truedrift),
  DIFFUSION = matrix(truediffusion), MANIFESTVAR = matrix(truemanifestvar),
  CINT = matrix(0), T0MEANS = matrix(0), T0VAR = matrix(1),
  MANIFESTMEANS = matrix(0))

traits <- stats::rnorm(nsub, 0, traitsd)
datalong <- do.call(rbind, lapply(seq_len(nsub), function(i) {
  subjectmodel <- generating
  subjectmodel$matrices$CINT <- matrix(traits[i])
  d <- ctGenerate(subjectmodel, n.subjects = 1, burnin = 20, dtmean = 1,
    Tpoints = tpoints, backend = 'r')
  d[, 'id'] <- i
  d
}))

estmodel <- ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
  LAMBDA = matrix(1), DRIFT = matrix(truedrift),
  DIFFUSION = matrix(truediffusion), MANIFESTVAR = matrix(truemanifestvar),
  CINT = matrix('cint1'), T0MEANS = matrix(0), T0VAR = matrix(1),
  MANIFESTMEANS = matrix(0))
estmodel$pars$indvarying[estmodel$pars$matrix == 'CINT'] <- TRUE

fit <- ctFit(datalong, estmodel, backend = 'julia', cores = 1, verbose = 0)

# A model whose DRIFT depends on the state, for the two routes' refusals and
# for the simulation route itself.
nonlinearmodel <- ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
  PARS = c('drift1'),
  LAMBDA = matrix(1), DRIFT = matrix('drift1 * exp(eta1 * 0.01)'),
  DIFFUSION = matrix(1), MANIFESTVAR = matrix(0.4), CINT = matrix(0),
  T0MEANS = matrix(0), T0VAR = matrix(1), MANIFESTMEANS = matrix(0))

# The same model with the random effect removed. Both routes can do this one,
# which is what lets them be checked against each other.
fixedmodel <- estmodel
fixedmodel$pars$indvarying <- FALSE
fixedfit <- ctFit(datalong, fixedmodel, backend = 'julia', cores = 1, verbose = 0)


test_that('the four components sum to the total and none is negative', {
  for (source in c('model', 'estimated')) {
    out <- ctVarianceDecomposition(fit, persons = source, npersons = 100)
    expect_equal(out$between + out$within.deterministic + out$within.stochastic +
        out$within.measurement, out$total)
    expect_true(all(out$between >= 0))
    expect_true(all(out$within.deterministic >= 0))
    expect_true(all(out$within.stochastic >= 0))
    expect_true(all(out$within.measurement >= 0))
    expect_equal(out$prop.between + out$prop.within, rep(1, nrow(out)))
  }
})

test_that('the stochastic and measurement parts match their closed forms', {
  out <- ctVarianceDecomposition(fit, persons = 'model', npersons = 60)
  matrices <- ctSummaryMatrices(fit)
  # The generating T0VAR is exactly the asymptotic variance, so the process
  # covariance is constant over the design and the stochastic part is that
  # value with no averaging error at all. Written out rather than read from the
  # fit: -DIFFUSIONcov / (2 DRIFT) for one process.
  asymptotic <- -as.numeric(matrices$DIFFUSIONcov) /
    (2 * as.numeric(matrices$DRIFT))
  expect_equal(out$within.stochastic[out$variable == 'eta1'], asymptotic)
  expect_equal(out$within.stochastic[out$variable == 'Y1'], asymptotic)
  expect_equal(out$within.measurement[out$variable == 'Y1'],
    as.numeric(matrices$MANIFESTcov))
  # A latent process has no measurement component at all.
  expect_equal(out$within.measurement[out$variable == 'eta1'], 0)
})

test_that('the deterministic part is the transient the model implies', {
  out <- ctVarianceDecomposition(fit, persons = 'model', npersons = 400)
  matrices <- ctSummaryMatrices(fit)
  drift <- as.numeric(matrices$DRIFT)
  # T0MEANS is fixed at zero and each person's asymptote is -CINT/DRIFT, so the
  # mean path is a(1 - exp(drift t)) and its variance over the design is
  # var_t(1 - exp(drift t)) * E[a^2]. Both factors are computed here rather
  # than taken from the function under test.
  times <- seq_len(tpoints) - 1
  shape <- 1 - exp(drift * times)
  shapevar <- mean((shape - mean(shape))^2)
  asymsd <- as.numeric(summary(fit)$popsd['cint1', 'mean']) / abs(drift)
  # E[a^2] is var + mean^2, and the population mean of cint1 is not exactly
  # zero even though it was generated that way.
  asymmean <- as.numeric(summary(fit)$popmeans['cint1', 'mean']) / abs(drift)
  expect_equal(out$within.deterministic[out$variable == 'eta1'],
    shapevar * (asymsd^2 + asymmean^2), tolerance = 0.15)
  # The between person variance is the same mean path averaged over time, so it
  # is var_t-free: mean(shape)^2 times the spread of the asymptotes.
  expect_equal(out$between[out$variable == 'eta1'],
    mean(shape)^2 * asymsd^2, tolerance = 0.15)
})

test_that('estimated persons are shrunk relative to drawn ones', {
  drawn <- ctVarianceDecomposition(fit, persons = 'model', npersons = 400)
  estimated <- ctVarianceDecomposition(fit, persons = 'estimated')
  expect_lt(estimated$between[estimated$variable == 'eta1'],
    drawn$between[drawn$variable == 'eta1'])
  # Shrinkage moves the persons, not the measurement model.
  expect_equal(estimated$within.measurement, drawn$within.measurement)
  expect_equal(attr(estimated, 'persons'), 'estimated')
  expect_equal(attr(drawn, 'npersons'), 400L)
})

test_that('a model with no random effects has no between person variance', {
  out <- ctVarianceDecomposition(fixedfit, persons = 'model', npersons = 20)
  # Every person has the same parameters and the same design, so the only
  # remaining source of between person variance is a design that differs, and
  # this one does not.
  expect_equal(out$between, rep(0, nrow(out)))
})

test_that('a state dependent model is refused by the moment route', {
  nonlinearfit <- ctFit(datalong[datalong[, 'id'] <= 6, ], nonlinearmodel,
    backend = 'julia', cores = 1, verbose = 0)
  expect_error(ctVarianceDecomposition(nonlinearfit, method = 'moment'),
    'depend on the latent state')
  # and the default picks the route that can answer.
  auto <- ctVarianceDecomposition(nonlinearfit, method = 'auto', npaths = 4)
  expect_equal(attr(auto, 'method'), 'simulation')
})


# The simulation route --------------------------------------------------------
#
# The check that matters is against the moment route on a model both can do.
# They share the estimand and nothing else: one propagates moments through
# expm(), the other draws trajectories in the engine and takes sample moments
# over them, so agreement is evidence about both rather than about neither.

test_that('simulation and moment agree on a model both can do', {
  set.seed(19)
  moment <- ctVarianceDecomposition(fixedfit, method = 'moment')
  simulated <- ctVarianceDecomposition(fixedfit, method = 'simulation',
    npaths = 400)
  expect_equal(attr(simulated, 'method'), 'simulation')
  expect_equal(simulated$within.measurement, moment$within.measurement)
  expect_equal(simulated$within.stochastic, moment$within.stochastic,
    tolerance = 0.05)
  expect_equal(simulated$total, moment$total, tolerance = 0.05)
  # Neither the process nor the design gives this model a between person
  # difference, and both routes have to say so.
  expect_equal(moment$between, rep(0, nrow(moment)))
  # The moment route gets exactly zero -- the model declares no levels, so
  # there is nothing to take a variance over. The simulation route estimates
  # the same zero from draws, so it is small rather than exact, and how small
  # depends on how many draws preceded it.
  expect_lt(max(simulated$between), 0.01 * min(simulated$total))
  expect_equal(simulated$between + simulated$within, simulated$total)
})

test_that('simulation and moment agree on a model with random effects', {
  # The harder half of the comparison above, and the one that exercises the
  # engine actually drawing the carrier states: with random effects there is a
  # between person term for the two routes to disagree about.
  set.seed(23)
  moment <- ctVarianceDecomposition(fit, method = 'moment', npersons = 600)
  simulated <- ctVarianceDecomposition(fit, method = 'simulation',
    npersons = 120, npaths = 30)
  expect_equal(simulated$within.measurement, moment$within.measurement)
  expect_equal(simulated$within.stochastic, moment$within.stochastic,
    tolerance = 0.05)
  expect_equal(simulated$between, moment$between, tolerance = 0.15)
  expect_equal(simulated$total, moment$total, tolerance = 0.05)
  # The point of the comparison: a between person variance that is actually
  # there. Before the engine drew the carrier states this came out at zero
  # while the moment route reported the real value, so what is checked is that
  # the two agree *and* that what they agree on is not zero. Scaled by the
  # moment route rather than given an absolute floor, which would depend on
  # whatever popsd this fixture's fit happens to land on -- 0.049 here, and a
  # floor of 0.1 guessed from another one is how this line first failed.
  expect_gt(min(moment$between), 1e-3)
  expect_gt(min(simulated$between), 0.5 * min(moment$between))
})

test_that('simulation refuses what it cannot answer', {
  expect_error(ctVarianceDecomposition(fixedfit, method = 'simulation',
    npaths = 1), 'npaths must be at least 2')
  expect_error(ctVarianceDecomposition(fixedfit, method = 'simulation',
    persons = 'estimated'), "persons='model'")
  expect_error(ctVarianceDecomposition(ctstantestfit, method = 'simulation'),
    "needs a backend='julia' fit")
})

test_that('the simulation route runs a state dependent model', {
  nonlinearfixed <- nonlinearmodel
  nonlinearfixed$pars$indvarying <- FALSE
  nonlinearfit <- ctFit(datalong, nonlinearfixed, backend = 'julia', cores = 1,
    verbose = 0)
  out <- ctVarianceDecomposition(nonlinearfit, npaths = 50)
  expect_equal(attr(out, 'method'), 'simulation')
  expect_equal(out$between + out$within, out$total)
  expect_true(all(out$within.stochastic > 0))
  expect_equal(out$within.measurement[out$type == 'latent'], 0)
  # The evaluation point is not reported because there is not one, and the
  # printed output has to say so rather than leave it implied.
  expect_output(print(out), 'drawn trajectories')
})


test_that('a binary indicator decomposes on both scales', {
  binarymodel <- ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
    LAMBDA = matrix(1), DRIFT = matrix(-0.5), DIFFUSION = matrix(1),
    MANIFESTVAR = matrix(0), CINT = matrix('cint1'), T0MEANS = matrix(0),
    T0VAR = matrix(1), MANIFESTMEANS = matrix(0), manifesttype = 1)
  binarymodel$pars$indvarying[binarymodel$pars$matrix == 'CINT'] <- TRUE
  binarydata <- datalong
  binarydata[, 'Y1'] <- as.numeric(binarydata[, 'Y1'] > 0)
  binaryfit <- ctFit(binarydata, binarymodel, backend = 'julia', cores = 1,
    verbose = 0)

  # The same seed for both, so that the two runs draw the same persons and any
  # difference between them is the scale rather than Monte Carlo error.
  set.seed(7)
  latent <- ctVarianceDecomposition(binaryfit, persons = 'model', npersons = 40,
    scale = 'latent')
  # The threshold form is a cumulative logit, so the response's own variance on
  # the latent scale is the logistic one.
  expect_equal(latent$within.measurement[latent$variable == 'Y1'], pi^2 / 3)

  set.seed(7)
  response <- ctVarianceDecomposition(binaryfit, persons = 'model',
    npersons = 40, scale = 'response')
  expect_equal(response$between + response$within, response$total)
  # A Bernoulli variable's variance cannot exceed 0.25, whatever the process
  # behind it does -- which the latent scale decomposition has no bound on.
  # Not asserted to the last digit: the between term is a sample variance over
  # the drawn persons and the time terms are averages over the design, so with
  # few persons the sum can sit a fraction of a percent above the bound it
  # converges to.
  expect_lt(response$total[response$variable == 'Y1'], 0.26)
  expect_gt(latent$total[latent$variable == 'Y1'], 3)
  # The scale argument moves the indicator and nothing else: the latent
  # processes behind it are the same quantity either way, and with the same
  # persons drawn they are the same numbers.
  numeric <- setdiff(names(latent), c('variable', 'type'))
  expect_equal(as.data.frame(latent)[latent$type == 'latent', numeric],
    as.data.frame(response)[response$type == 'latent', numeric])
})


test_that('a time dependent predictor shows up as deterministic within variance', {
  tdmodel <- ctModel(type = 'ct', n.latent = 1, n.manifest = 1, n.TDpred = 1,
    LAMBDA = matrix(1), DRIFT = matrix(-0.5), DIFFUSION = matrix(1),
    MANIFESTVAR = matrix(0.4), CINT = matrix(0), T0MEANS = matrix(0),
    T0VAR = matrix(1), MANIFESTMEANS = matrix(0),
    TDPREDEFFECT = matrix('tdeffect'))
  tddata <- cbind(datalong, TD1 = 0)
  tddata[tddata[, 'time'] %in% c(5, 6, 7), 'TD1'] <- 2
  tdfit <- ctFit(tddata, tdmodel, backend = 'julia', cores = 1, verbose = 0)
  out <- ctVarianceDecomposition(tdfit, persons = 'model', npersons = 20)
  expect_gt(out$within.deterministic[out$variable == 'eta1'], 0)
  expect_equal(out$between + out$within, out$total)
})


test_that('the two backends decompose the same model the same way', {
  # Only when CTSEM_TEST_STAN asks for it, as the rest of this suite does: the
  # stan fit is the expensive half. What it checks is the one thing the two
  # paths do differently -- stan reads the subject matrices the filter saved and
  # falls back to pop_ for the rest, julia materialises each person at its
  # carrier vector -- so agreement here is about those two routes rather than
  # about the decomposition, which is shared.
  skip_if(!identical(test_backends(), c('julia', 'stan')),
    'set CTSEM_TEST_STAN to compare the backends')
  fits <- fit_backends(datalong = datalong, model = estmodel, cores = 1,
    verbose = 0)
  julia <- ctVarianceDecomposition(fits$julia, persons = 'estimated')
  stan <- suppressMessages(ctVarianceDecomposition(fits$stan, persons = 'estimated'))
  parts <- c('between', 'within.deterministic', 'within.stochastic',
    'within.measurement', 'total')
  for (part in parts) {
    expect_equal(stan[[part]], julia[[part]], tolerance = 0.05,
      label = paste0('stan ', part))
  }
})

# Resolved over time ----------------------------------------------------------
#
# The aggregate is a window average, and for a process that has not settled the
# split it averages is not the split at any particular occasion. `times=`
# evaluates every person on one common grid instead.

test_that('times= resolves the decomposition over the window', {
  grid <- c(0, 1, 2, 5, 10, 19)
  out <- ctVarianceDecomposition(fit, times = grid, npersons = 300)

  expect_true('time' %in% names(out))
  expect_equal(nrow(out), length(grid) * 2L)     # two variables
  expect_equal(sort(unique(out$time)), grid)
  # There is no deterministic term at an instant: it is the variance of the
  # mean path *over* time.
  expect_false('within.deterministic' %in% names(out))
  expect_equal(out$between + out$within, out$total)

  # Each row must be about the variable it names. A transpose in the flattening
  # gives a table where every row still sums and none of them is, which is how
  # this was first written.
  expect_true(all(out$within.measurement[out$type == 'latent'] == 0))
  expect_true(all(out$within.measurement[out$type == 'manifest'] > 0))

  # T0MEANS is fixed, so at the first occasion every person starts in the same
  # place and there is no between person variance at all; it accumulates as the
  # trajectories separate toward their own asymptotes.
  # Stated without a magnitude. The claim is that the split changes across the
  # window, and how far it gets depends on this fixture's popsd and occasion
  # count -- a floor read off some other fixture is how the first three
  # versions of this file failed.
  first <- out[out$time == min(grid) & out$variable == 'eta1', ]
  last <- out[out$time == max(grid) & out$variable == 'eta1', ]
  expect_equal(first$between, 0)
  expect_gt(last$between, 0)
  expect_gt(last$prop.between, first$prop.between)
  # and it is monotone here, which is what a relaxation looks like.
  series <- out$between[out$variable == 'eta1'][order(unique(out$time))]
  expect_true(all(diff(series) >= 0))
})

test_that('the two modes agree about the part they share', {
  # The stochastic term is a mean over occasions either way, so evaluating it
  # on the design's own times and averaging must give the window average. The
  # between term is not comparable -- one is the variance of a time average,
  # the other a time series of variances -- which is the whole reason both
  # exist.
  set.seed(31)
  aggregate <- ctVarianceDecomposition(fit, npersons = 300)
  set.seed(31)
  resolved <- ctVarianceDecomposition(fit, times = 'asdata', npersons = 300)
  bytime <- tapply(resolved$within.stochastic[resolved$variable == 'eta1'],
    resolved$time[resolved$variable == 'eta1'], mean)
  expect_equal(mean(bytime),
    aggregate$within.stochastic[aggregate$variable == 'eta1'],
    tolerance = 1e-8)
  expect_equal(unique(resolved$time), sort(unique(datalong[, 'time'])))
})

test_that('times= refuses what it has no grid for', {
  expect_error(ctVarianceDecomposition(fit, times = c(0, 1),
    method = 'simulation'), 'no common grid')
  expect_error(ctVarianceDecomposition(fit, times = numeric()), 'empty')
})


# Designs and levels ----------------------------------------------------------

test_that('the decomposition holds together on any design shape', {
  # Balanced, irregular, ragged and missing, on one fit each. What is checked is
  # the identity and that nothing comes back NA -- the values differ between
  # them because the fits do, and comparing them would be comparing fits.
  shapes <- list(
    balanced = lapply(1:10, function(i) 0:9),
    irregular = lapply(1:10, function(i) sort(cumsum(c(0, stats::rexp(9, 1))))),
    ragged = lapply(1:10, function(i) 0:(2 + i)))
  set.seed(53)
  traits <- stats::rnorm(10, 0, 0.35)
  generating <- estmodel
  generating$pars$indvarying <- FALSE

  for (shape in names(shapes)) {
    times <- shapes[[shape]]
    dat <- do.call(rbind, lapply(seq_along(traits), function(i) {
      m <- generating
      m$matrices$CINT <- matrix(traits[i])
      d <- suppressMessages(ctGenerate(m, n.subjects = 1, burnin = 15,
        Tpoints = length(times[[i]]), backend = 'r'))
      d <- as.data.frame(d)
      d$time <- times[[i]]
      d$id <- i
      d
    }))
    if (identical(shape, 'balanced')) {
      dat$Y1[seq(1, nrow(dat), by = 3)] <- NA   # and a third of it missing
    }
    f <- suppressWarnings(suppressMessages(ctFit(dat, estmodel,
      backend = 'julia', cores = 1, verbose = 0)))
    out <- suppressMessages(ctVarianceDecomposition(f, npersons = 60))
    expect_equal(out$between + out$within, out$total, label = shape)
    expect_false(anyNA(out$total))
    expect_true(all(out$within.stochastic > 0))
  }
})

test_that('the two representations of individual differences agree', {
  # The test that would have caught the Laplace gap in one line, and the reason
  # `fit_intoverpop()` exists: the same model, the same data, the same
  # individual differences, described in two places that share no field. Read
  # only one of them and the other returns a between person variance of zero --
  # which is what a model with no individual differences correctly returns, so
  # nothing about the number looks wrong.
  #
  # Compared at persons='estimated' on both sides, since that is the only route
  # a Laplace fit has: both are then the fitted modes and the same quantity.
  fits <- fit_intoverpop(datalong = datalong, model = estmodel, cores = 1,
    verbose = 0)
  augmented <- suppressMessages(
    ctVarianceDecomposition(fits$augmented, persons = 'estimated'))
  laplace <- suppressMessages(
    ctVarianceDecomposition(fits$laplace, persons = 'estimated'))

  expect_gt(min(augmented$between), 0)
  expect_gt(min(laplace$between), 0)
  # The same body over the same person structure, so at persons='estimated'
  # -- where neither draws anything -- the components differ only as far as
  # the two fits' estimates do. Those are separate optimisations, each stopped
  # by its own rule (the augmented one finishes with Newton steps, the laplace
  # one does not), and agree to about 5e-6 relative. It was the between term
  # differing by two thirds that started this, so 1e-4 still says 'the same'.
  for (part in c('between', 'within.deterministic', 'within.stochastic',
    'within.measurement', 'total')) {
    expect_equal(laplace[[part]], augmented[[part]], tolerance = 1e-4,
      label = paste('laplace', part))
  }
})

test_that('drawing persons agrees across representations too', {
  # The same comparison with both sides drawing from their own fitted
  # population covariance rather than reading modes. Looser, because the two
  # draw independently, but it is the route the default takes.
  fits <- fit_intoverpop(datalong = datalong, model = estmodel, cores = 1,
    verbose = 0)
  drawn <- lapply(fits, function(f) {
    set.seed(29)
    suppressMessages(ctVarianceDecomposition(f, persons = 'model',
      npersons = 800))
  })
  for (part in c('between', 'within.stochastic', 'within.measurement',
    'total')) {
    expect_equal(drawn$laplace[[part]], drawn$augmented[[part]],
      tolerance = 0.12, label = paste('drawn', part))
  }
  # And drawing is not shrunk, so it exceeds what the modes give.
  estimated <- suppressMessages(
    ctVarianceDecomposition(fits$laplace, persons = 'estimated'))
  expect_gt(min(drawn$laplace$between), min(estimated$between))
})

test_that('the random effect structure reads the same from either representation', {
  # The accessor the decomposition now asks instead of looking for carrier
  # states, over both representations and over a model that has none.
  fits <- fit_intoverpop(datalong = datalong, model = estmodel, cores = 1,
    verbose = 0)
  for (rep in names(fits)) {
    levels <- ctsem:::.ctFitRandomEffectLevels(fits[[rep]])
    expect_length(levels, 1L)
    expect_equal(levels[[1L]]$params, 'cint1', label = rep)
    expect_equal(levels[[1L]]$nunits, nsub, label = rep)
    expect_equal(levels[[1L]]$units, seq_len(nsub), label = rep)
    expect_true(ctsem:::.ctFitHasRandomEffects(fits[[rep]]), label = rep)
  }
  # No individual differences is an empty answer, not an unreadable one -- the
  # distinction the Laplace failure turned on.
  expect_length(ctsem:::.ctFitRandomEffectLevels(fixedfit), 0L)
  expect_false(ctsem:::.ctFitHasRandomEffects(fixedfit))
  # And a stan fit answers the same question.
  expect_true(ctsem:::.ctFitHasRandomEffects(ctstantestfit))
})

test_that('a level above the subject gets its own between column', {
  # Only intoverpop='laplace' can integrate out a level above the subject, and
  # a Laplace fit has no carrier states -- which is exactly why this needs its
  # own route: the carrier one would report a between person variance of zero
  # for every multilevel model and say nothing.
  set.seed(61)
  generating <- estmodel
  generating$pars$indvarying <- FALSE
  ngroup <- 4L
  groupeffect <- stats::rnorm(ngroup, 0, 0.5)
  rows <- list()
  unit <- 0L
  for (g in seq_len(ngroup)) for (p in 1:4) {
    unit <- unit + 1L
    m <- generating
    m$matrices$CINT <- matrix(groupeffect[g] + stats::rnorm(1, 0, 0.3))
    d <- suppressMessages(ctGenerate(m, n.subjects = 1, burnin = 15,
      Tpoints = 10, backend = 'r'))
    d <- as.data.frame(d)
    d$id <- unit
    d$study <- g
    rows[[unit]] <- d
  }
  dat <- do.call(rbind, rows)

  model <- ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
    id = c('id', 'study'),
    LAMBDA = matrix(1), DRIFT = matrix(-0.5), DIFFUSION = matrix(1),
    MANIFESTVAR = matrix(0.4), CINT = matrix('cint1'), T0MEANS = matrix(0),
    T0VAR = matrix(1), MANIFESTMEANS = matrix(0))
  model$pars$indvarying <- model$pars$matrix == 'CINT'
  model$pars$indvarying_study <- model$pars$matrix == 'CINT'
  fit <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = 'julia',
    intoverpop = 'laplace', cores = 1, verbose = 0)))

  out <- suppressMessages(ctVarianceDecomposition(fit))
  expect_true(all(c('between.id', 'between.study') %in% names(out)))
  # The levels are orthogonal by construction, so they add up to `between`.
  expect_equal(out$between.id + out$between.study, out$between)
  expect_equal(out$between + out$within, out$total)
  expect_true(all(out$between.study >= 0))
  expect_true(all(out$between.id >= 0))
  # A between person variance that is actually there: the carrier route this
  # replaces returned zero for the same fit.
  expect_gt(min(out$between), 0)
  # The drawn route works here too, and is not shrunk: a level's modes are
  # pulled toward the level outside them, and drawing from that level's own
  # fitted covariance is not. The subject level is where it shows, since a
  # study aggregates several subjects and is better determined than any of
  # them -- shrinkage follows the information per unit, not the unit count.
  drawn <- suppressMessages(ctVarianceDecomposition(fit, persons = 'model',
    npersons = 400))
  expect_true(all(c('between.id', 'between.study') %in% names(drawn)))
  expect_equal(drawn$between.id + drawn$between.study, drawn$between)
  expect_equal(drawn$between + drawn$within, drawn$total)
  expect_gt(min(drawn$between.id), min(out$between.id))
})

test_that('a stan fit works through the estimated route and refuses the drawn one', {
  skip_on_cran()
  skip_on_32bit()
  out <- expect_message(ctVarianceDecomposition(ctstantestfit),
    "persons='estimated'")
  expect_equal(out$between + out$within, out$total)
  expect_true(all(out$total > 0))
  expect_error(ctVarianceDecomposition(ctstantestfit, persons = 'model'),
    "needs a backend='julia' fit")
})
