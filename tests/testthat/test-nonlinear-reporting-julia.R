library(ctsem)
library(testthat)

context('nonlinearreportingjulia')

# One nonlinear julia fit, reused by everything below: fitting is the expensive
# part and the questions here are all about what gets *reported* from it.
nlfit <- local({
  skip_without_julia()
  set.seed(1)
  generating <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
    LAMBDA = diag(2), DRIFT = matrix(c(-.4, .1, 0, -.3), 2, 2),
    CINT = matrix(c(.2, .1), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(.2, 2), DIFFUSION = matrix(c(.5, 0, 0, .4), 2, 2)))
  datalong <- as.data.frame(suppressMessages(ctGenerate(generating, n.subjects = 25,
    burnin = 5, dtmean = 1, logdtsd = .1, wide = FALSE, Tpoints = 12)))

  model <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
    LAMBDA = diag(2),
    # eta1's own decay depends on where eta2 is.
    PARS = c('dr11|-log1p_exp(param)'),
    DRIFT = matrix(c('dr11 * (1 + 0.2 * eta2)', 'd21', 0, 'd22'), 2, 2),
    CINT = matrix(c('c1', 'c2'), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(.2, 2), DIFFUSION = matrix(c('df1', 0, 0, 'df2'), 2, 2)))
  model$pars$indvarying <- FALSE
  ctFit(datalong, model, backend = 'julia', cores = 1, verbose = 0)
})

test_that('the state-dependent cell is found, and only it', {
  cells <- ctsem:::.ctFitContextDependentCells(nlfit)
  expect_true(ctModelIsNonlinear(nlfit))
  expect_true(all(cells$kind == 'state'))
  expect_true('DRIFT' %in% cells$matrix)
  expect_true(any(cells$matrix == 'DRIFT' & cells$row == 1 & cells$col == 1))
  # DRIFT[2,2] is a plain parameter and must not be swept in.
  expect_false(any(cells$matrix == 'DRIFT' & cells$row == 2 & cells$col == 2))
})

test_that('the evaluation point changes the dependent cell and nothing else', {
  t0 <- suppressMessages(ctSummaryMatrices(nlfit))
  mean <- suppressMessages(ctSummaryMatrices(nlfit, state = 'mean'))
  asymptotic <- suppressMessages(ctSummaryMatrices(nlfit, state = 'asymptotic'))
  explicit <- suppressMessages(ctSummaryMatrices(nlfit, state = c(0, 5)))

  expect_false(isTRUE(all.equal(t0$DRIFT[1, 1], mean$DRIFT[1, 1])))
  expect_false(isTRUE(all.equal(t0$DRIFT[1, 1], explicit$DRIFT[1, 1])))
  # Cells that do not depend on the state must be untouched by the choice.
  expect_equal(t0$DRIFT[2, 2], mean$DRIFT[2, 2])
  expect_equal(t0$DRIFT[2, 2], asymptotic$DRIFT[2, 2])
  expect_equal(t0$DIFFUSIONcov, mean$DIFFUSIONcov)

  # DRIFT[1,1] is dr11 * (1 + 0.2 * eta2), so doubling the multiplier by
  # setting eta2 = 5 must double the cell relative to eta2 = 0.
  atzero <- suppressMessages(ctSummaryMatrices(nlfit, state = c(0, 0)))
  expect_equal(explicit$DRIFT[1, 1], 2 * atzero$DRIFT[1, 1], tolerance = 1e-6)
})

test_that('the summary names the point, and does not name Jacobian blocks', {
  note <- summary(nlfit)$parmatNote
  expect_match(note, 'DRIFT')
  expect_match(note, 'T0MEANS state')
  expect_false(grepl('JAx', note, fixed = TRUE))
})

test_that('a supplied state is validated and padded', {
  expect_error(suppressMessages(ctSummaryMatrices(nlfit, state = c(1, 2, 3, 4, 5))),
    'must have 2 entries')
  expect_length(ctsem:::.ctResolveState(nlfit, 'mean')$state,
    length(ctsem:::.ctContextBaseState(nlfit)))
})

test_that("method='simulate' reduces exactly to the linearised answer when linear", {
  linear <- local({
    set.seed(1)
    generating <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
      manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
      LAMBDA = diag(2), DRIFT = matrix(c(-.4, .1, 0, -.3), 2, 2),
      CINT = matrix(c(.2, .1), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
      MANIFESTVAR = diag(.2, 2), DIFFUSION = matrix(c(.5, 0, 0, .4), 2, 2)))
    datalong <- as.data.frame(suppressMessages(ctGenerate(generating,
      n.subjects = 25, burnin = 5, dtmean = 1, logdtsd = .1, wide = FALSE,
      Tpoints = 12)))
    model <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
      manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
      LAMBDA = diag(2), DRIFT = matrix(c('d11', 'd21', 0, 'd22'), 2, 2),
      CINT = matrix(c('c1', 'c2'), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
      MANIFESTVAR = diag(.2, 2), DIFFUSION = matrix(c('df1', 0, 0, 'df2'), 2, 2)))
    model$pars$indvarying <- FALSE
    ctFit(datalong, model, backend = 'julia', cores = 1, verbose = 0)
  })
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
  simulated <- suppressMessages(ctDiscretePars(nlfit, times = times,
    method = 'simulate', nsamples = 3))
  linearised <- suppressMessages(ctDiscretePars(nlfit, times = times,
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
