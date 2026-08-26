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
