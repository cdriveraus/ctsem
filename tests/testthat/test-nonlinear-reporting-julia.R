library(ctsem)
library(testthat)

context('nonlinearreportingjulia')

# Fitting is by far the expensive part of this file and none of the questions
# here are about fitting, so each model is built once, on first use, and shared.
# Two fits cover everything: one state dependent, one linear with correlated
# diffusion and unequal process scales -- the linear one has to be correlated,
# because a diagonal diffusion never exercises the companion-shock path at all
# and would let a wrong companion matrix pass every reduction check.
.ctTestFit <- local({
  store <- list()
  function(name, builder) {
    if (is.null(store[[name]])) store[[name]] <<- builder()
    store[[name]]
  }
})

.ctTestData <- function(seed, diffusion, nsubjects = 25, Tpoints = 12) {
  set.seed(seed)
  generating <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
    LAMBDA = diag(2), DRIFT = matrix(c(-.4, .1, 0, -.3), 2, 2),
    CINT = matrix(c(.2, .1), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(.2, 2), DIFFUSION = diffusion))
  as.data.frame(suppressMessages(ctGenerate(generating, n.subjects = nsubjects,
    burnin = 5, dtmean = 1, logdtsd = .1, wide = FALSE, Tpoints = Tpoints)))
}

# State dependent: eta1's own decay depends on where eta2 is.
nonlinearFit <- function() .ctTestFit('nonlinear', function() {
  skip_without_julia()
  datalong <- .ctTestData(1, matrix(c(.5, 0, 0, .4), 2, 2))
  model <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
    LAMBDA = diag(2), PARS = c('dr11|-log1p_exp(param)'),
    DRIFT = matrix(c('dr11 * (1 + 0.2 * eta2)', 'd21', 0, 'd22'), 2, 2),
    CINT = matrix(c('c1', 'c2'), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(.2, 2), DIFFUSION = matrix(c('df1', 0, 0, 'df2'), 2, 2)))
  model$pars$indvarying <- FALSE
  ctFit(datalong, model, backend = 'julia', cores = 1, verbose = 0)
})

# Linear, with correlated diffusion and processes on different scales.
linearFit <- function() .ctTestFit('linear', function() {
  skip_without_julia()
  datalong <- .ctTestData(5, matrix(c(1.0, 0.8, 0, 1.2), 2, 2), nsubjects = 30)
  model <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
    LAMBDA = diag(2), DRIFT = matrix(c('d11', 'd21', 0, 'd22'), 2, 2),
    CINT = matrix(c('c1', 'c2'), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(.2, 2),
    DIFFUSION = matrix(c('df11', 'df21', 0, 'df22'), 2, 2)))
  model$pars$indvarying <- FALSE
  ctFit(datalong, model, backend = 'julia', cores = 1, verbose = 0)
})

# State dependent, with a time independent predictor, for the covariate panel.
tipredFit <- function() .ctTestFit('tipred', function() {
  skip_without_julia()
  datalong <- .ctTestData(2, matrix(c(.5, 0, 0, .4), 2, 2), nsubjects = 30, Tpoints = 10)
  set.seed(3)
  values <- rnorm(length(unique(datalong$id)))
  datalong$TI1 <- values[match(datalong$id, unique(datalong$id))]
  model <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    n.TIpred = 1, manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
    TIpredNames = 'TI1', LAMBDA = diag(2), PARS = c('dr11|-log1p_exp(param)'),
    DRIFT = matrix(c('dr11 * (1 + 0.2 * eta2)', 'd21', 0, 'd22'), 2, 2),
    CINT = matrix(c('c1', 'c2'), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(.2, 2), DIFFUSION = matrix(c('df1', 0, 0, 'df2'), 2, 2)))
  model$pars$indvarying <- FALSE
  ctFit(datalong, model, backend = 'julia', cores = 1, verbose = 0)
})


test_that('the state-dependent cell is found, and only it', {
  cells <- ctsem:::.ctFitContextDependentCells(nonlinearFit())
  expect_true(ctModelIsNonlinear(nonlinearFit()))
  expect_true(all(cells$kind == 'state'))
  expect_true('DRIFT' %in% cells$matrix)
  expect_true(any(cells$matrix == 'DRIFT' & cells$row == 1 & cells$col == 1))
  # DRIFT[2,2] is a plain parameter and must not be swept in.
  expect_false(any(cells$matrix == 'DRIFT' & cells$row == 2 & cells$col == 2))
})

test_that('the evaluation point changes the dependent cell and nothing else', {
  t0 <- suppressMessages(ctSummaryMatrices(nonlinearFit()))
  mean <- suppressMessages(ctSummaryMatrices(nonlinearFit(), state = 'mean'))
  asymptotic <- suppressMessages(ctSummaryMatrices(nonlinearFit(), state = 'asymptotic'))
  explicit <- suppressMessages(ctSummaryMatrices(nonlinearFit(), state = c(0, 5)))

  expect_false(isTRUE(all.equal(t0$DRIFT[1, 1], mean$DRIFT[1, 1])))
  expect_false(isTRUE(all.equal(t0$DRIFT[1, 1], explicit$DRIFT[1, 1])))
  # Cells that do not depend on the state must be untouched by the choice.
  expect_equal(t0$DRIFT[2, 2], mean$DRIFT[2, 2])
  expect_equal(t0$DRIFT[2, 2], asymptotic$DRIFT[2, 2])
  expect_equal(t0$DIFFUSIONcov, mean$DIFFUSIONcov)

  # DRIFT[1,1] is dr11 * (1 + 0.2 * eta2), so doubling the multiplier by
  # setting eta2 = 5 must double the cell relative to eta2 = 0.
  atzero <- suppressMessages(ctSummaryMatrices(nonlinearFit(), state = c(0, 0)))
  expect_equal(explicit$DRIFT[1, 1], 2 * atzero$DRIFT[1, 1], tolerance = 1e-6)
})

# The julia half of ctExtract()'s state= argument; the stan half (which refuses
# it) is in test-context-dependence.R. Before this, the shorthands reached the
# engine as a bare string and died there with a Julia MethodError, an explicit
# vector worked but came back with nothing saying where it had been evaluated,
# and the stan method dropped the argument entirely.
test_that('ctExtract honours state= on julia and records the point', {
  fit <- nonlinearFit()
  default <- suppressMessages(ctExtract(fit))
  mean <- suppressMessages(ctExtract(fit, state = 'mean'))
  asymptotic <- suppressMessages(ctExtract(fit, state = 'asymptotic'))
  atzero <- suppressMessages(ctExtract(fit, state = c(0, 0)))
  atfive <- suppressMessages(ctExtract(fit, state = c(0, 5)))

  # Same words as ctSummaryMatrices() uses for the same points.
  expect_identical(attr(default, 'evaluatedAt'), ctsem:::.ctContextPopLabel)
  expect_identical(attr(mean, 'evaluatedAt'), 'the mean smoothed latent state')
  expect_identical(attr(asymptotic, 'evaluatedAt'),
    "the system's asymptotic (fixed point) state")
  expect_identical(attr(atfive, 'evaluatedAt'), 'the supplied state')

  # The state dependent cell moves and the others do not -- the same claim
  # ctSummaryMatrices() makes, checked on the arrays it is computed from.
  expect_false(isTRUE(all.equal(default$pop_DRIFT[, 1, 1], mean$pop_DRIFT[, 1, 1])))
  expect_equal(default$pop_DRIFT[, 2, 2], mean$pop_DRIFT[, 2, 2])
  expect_equal(default$pop_DIFFUSION, mean$pop_DIFFUSION)
  # DRIFT[1,1] is dr11 * (1 + 0.2 * eta2): eta2 = 5 doubles it against eta2 = 0.
  expect_equal(atfive$pop_DRIFT[, 1, 1], 2 * atzero$pop_DRIFT[, 1, 1], tolerance = 1e-6)

  # Recorded, never messaged: ctExtract() is called repeatedly by the summaries.
  expect_no_message(ctExtract(fit))
  expect_no_message(ctExtract(fit, state = 'mean'))
})

test_that('the summary names the point, and does not name Jacobian blocks', {
  note <- summary(nonlinearFit())$parmatNote
  expect_match(note, 'DRIFT')
  expect_match(note, 'T0MEANS state')
  expect_false(grepl('JAx', note, fixed = TRUE))
})

test_that('a supplied state is validated and padded', {
  expect_error(suppressMessages(ctSummaryMatrices(nonlinearFit(), state = c(1, 2, 3, 4, 5))),
    'must have 2 entries')
  expect_length(ctsem:::.ctResolveState(nonlinearFit(), 'mean')$state,
    length(ctsem:::.ctContextBaseState(nonlinearFit())))
})

test_that("method='simulate' reduces exactly to the linearised answer when linear", {
  linear <- linearFit()
  expect_false(ctModelIsNonlinear(linear))

  times <- c(0, .5, 1, 2, 4)
  simulated <- suppressMessages(ctDiscretePars(linear, times = times,
    method = 'simulate', nsamples = 3))
  linearised <- suppressMessages(ctDiscretePars(linear, times = times, nsamples = 3))
  # Integrating a linear system and exponentiating its drift are the same
  # calculation; if they are not, the integrator is wrong.
  expect_equal(apply(simulated, c(3, 4, 5), median),
    apply(linearised, c(3, 4, 5), median), tolerance = 1e-10)

  # And the integrator against the closed form directly, at the point estimate.
  state <- ctsem:::.ctContextBaseState(linear)
  response <- ctsem:::.ctNonlinearImpulseResponse(linear, state, times)
  drift <- suppressMessages(ctBackendParMatrices(linear, trim = FALSE))$DRIFT
  for (k in seq_along(times)) {
    expect_equal(unname(response[k, , ]),
      unname(as.matrix(Matrix::expm(drift * times[k]))), tolerance = 1e-10)
  }
})

test_that("method='simulate' finds an effect the linearisation reports as zero", {
  times <- c(0, 1, 4)
  simulated <- suppressMessages(ctDiscretePars(nonlinearFit(), times = times,
    method = 'simulate', nsamples = 3))
  linearised <- suppressMessages(ctDiscretePars(nonlinearFit(), times = times,
    state = 'asymptotic', nsamples = 3))

  # The impulse response starts at the identity by construction.
  expect_equal(unname(apply(simulated, c(3, 4, 5), median)[1, , ]), diag(2),
    tolerance = 1e-8)

  # DRIFT[1,2] is fixed at 0, so the linearised panel says eta2 can never affect
  # eta1. It can: eta1's decay *rate* depends on eta2, which no frozen DRIFT
  # can express. This is the whole reason the method exists.
  expect_equal(apply(linearised, c(3, 4, 5), median)[3, 1, 2], 0)
  expect_true(abs(apply(simulated, c(3, 4, 5), median)[3, 1, 2]) > 1e-3)
})

test_that("method='simulate' is refused where it cannot be honoured", {
  expect_error(ctDiscretePars(ctstantestfit, times = 1, method = 'simulate'),
    "backend='julia'")
})

test_that('the covariate dynamics panel is built from each level\'s own state', {
  tipfit <- tipredFit()

  panel <- suppressMessages(ctsem:::.ctPredictTIPDynamics(tipfit, tipredIndex = 1,
    values = c(-1, 0, 1), times = c(0, 1, 2, 4), ntipred = 1, nsamples = 3,
    latentNames = c('eta1', 'eta2')))
  median <- apply(panel, c(2, 3, 4, 5), stats::median)

  # The reading the linear panel has, preserved: a one unit impulse starts at
  # one and decays to zero for a stable process -- at every covariate level,
  # even though each level is evaluated at a different state.
  for (level in 1:3) expect_equal(unname(median[level, 1, , ]), diag(2),
    tolerance = 1e-8)
  expect_true(all(abs(median[, 4, 1, 1]) < abs(median[, 1, 1, 1])))

  # And the levels differ, which is what the panel is for.
  expect_true(diff(range(median[, 3, 1, 1])) > 1e-6)
})

test_that('every observational and standardise combination reduces when linear', {
  # Correlated diffusion and unequal process scales: the case that exposes a
  # correlation being used where a regression coefficient belongs. Both earlier
  # test models have diagonal diffusion and never exercise the companion path.
  correlated <- linearFit()

  mats <- suppressMessages(ctSummaryMatrices(correlated))
  expect_gt(abs(stats::cov2cor(mats$DIFFUSIONcov)[2, 1]), .2)
  scales <- sqrt(diag(mats$asymDIFFUSIONcov))
  expect_gt(max(scales) / min(scales), 1.3)

  times <- c(0, .5, 1, 2, 4)
  # Every interpretation, both scalings. Simulating a linear system and
  # exponentiating its drift are the same calculation whichever companion
  # matrix is applied, so any of these failing means the two paths disagree
  # about what was asked for.
  for (observational in ctsem:::.ctCompanionTypes) for (standardise in c(FALSE, TRUE)) {
    simulated <- suppressMessages(ctDiscretePars(correlated, times = times,
      method = 'simulate', nsamples = 3, observational = observational,
      standardise = standardise))
    linearised <- suppressMessages(ctDiscretePars(correlated, times = times,
      nsamples = 3, observational = observational, standardise = standardise))
    expect_equal(apply(simulated, c(3, 4, 5), stats::median),
      apply(linearised, c(3, 4, 5), stats::median), tolerance = 1e-9)
  }
})

test_that('the phase portrait recovers the fixed point and the curved nullcline', {
  portrait <- suppressMessages(ctPhasePortrait(nonlinearFit(), latents = c('eta1', 'eta2'),
    gridsize = 11, extent = 'sd', plot = FALSE))
  expect_true(all(is.finite(portrait$field$dx)))

  # The marked point is where the field vanishes, to solver tolerance.
  field <- ctsem:::.ctFieldFunction(nonlinearFit())
  state <- ctsem:::.ctResolveState(nonlinearFit(), 'asymptotic')$state
  expect_lt(max(abs(field(state))), 1e-6)

  # DRIFT[1,1] depends on eta2 while DRIFT[2,] does not depend on anything, so
  # eta1's nullcline must bend and eta2's must not. Straightness is measured as
  # the residual of a straight line fit through the contour.
  straightness <- vapply(split(portrait$nullclines, portrait$nullclines$process),
    function(piece) {
      if (nrow(piece) < 4) return(NA_real_)
      max(abs(stats::residuals(stats::lm(y ~ x, data = piece)))) /
        max(diff(range(piece$y)), 1e-8)
    }, numeric(1))
  expect_lt(straightness[['eta2']], 1e-6)
  expect_gt(straightness[['eta1']], straightness[['eta2']])
})

test_that('the state dependence plot recovers the specified relationship', {
  # DRIFT[1,1] was specified as dr11 * (1 + 0.2 * eta2), so its value must be
  # exactly linear in eta2 with slope 0.2 * dr11. Recovering that from the
  # fitted engine is a check on the whole materialise-at-a-state path.
  values <- suppressMessages(ctStateDependencePlot(nonlinearFit(), along = 'eta2',
    gridsize = 11, nsamples = 5, plot = FALSE))
  expect_equal(unique(values$cell), 'DRIFT[1,1]')

  straight <- stats::lm(middle ~ along, data = values)
  # suppressWarnings: R objects to an essentially perfect fit, which is the
  # property being tested.
  expect_gt(suppressWarnings(summary(straight))$r.squared, 1 - 1e-9)
  intercept <- unname(stats::coef(straight)[1])
  slope <- unname(stats::coef(straight)[2])
  expect_equal(slope, 0.2 * intercept, tolerance = 1e-6)
})

test_that('the state dependence plot declines a linear model rather than drawing nothing', {
  linear <- linearFit()
  expect_error(ctStateDependencePlot(linear), 'nothing to plot')
})
