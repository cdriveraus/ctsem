library(ctsem)
library(testthat)

context('contextdependence')

# The expression parser is the whole of the detection logic, so it is tested
# directly rather than only through a fit: a fit is expensive and would hide
# which half broke.

test_that('expression indices are extracted in every written form', {
  expect_equal(ctsem:::.ctExpressionIndices('2*state[3] + state[11]', 'state'), c(3L, 11L))
  expect_equal(ctsem:::.ctExpressionIndices('a*tdpreds[rowi, 2]', 'tdpreds'), 2L)
  expect_equal(ctsem:::.ctExpressionIndices('a*ctx.tdpreds[2]', 'tdpreds'), 2L)
  expect_length(ctsem:::.ctExpressionIndices('plainparam', 'state'), 0)
})

test_that('a state reference is classified by whether it is a carrier', {
  nlatent <- 2
  expect_equal(ctsem:::.ctExpressionDependence('-log1p(exp(param))*state[2]', nlatent), 'state')
  expect_equal(ctsem:::.ctExpressionDependence('param*state[5]', nlatent), 'carrier')
  expect_equal(sort(ctsem:::.ctExpressionDependence('state[1]*state[7]', nlatent)),
    c('carrier', 'state'))
  expect_equal(ctsem:::.ctExpressionDependence('b*tdpreds[rowi, 1]', nlatent), 'tdpred')
  expect_length(ctsem:::.ctExpressionDependence('drift11', nlatent), 0)
  expect_length(ctsem:::.ctExpressionDependence(NA_character_, nlatent), 0)
})

test_that('dependence propagates through PARS references', {
  fit <- structure(list(
    standata = list(nlatent = 2L),
    setup = list(matsetup = data.frame(
      parname = c('PARS[1,1] * param', 'log1p(exp(state[1]))', 'd22', 'param*state[9]'),
      row = c(1L, 1L, 2L, 1L), col = c(1L, 1L, 2L, 1L),
      matrix = c(3L, 10L, 3L, 7L)))),
    class = c('ctStanFit', 'ctFit'))

  cells <- ctsem:::.ctFitContextDependentCells(fit)
  expect_equal(nrow(cells), 3)
  # DRIFT[1,1] references PARS[1,1], which references state[1]: transitive.
  expect_true('state' %in% cells$kind[cells$matrix == 'DRIFT' & cells$row == 1])
  expect_true('carrier' %in% cells$kind[cells$matrix == 'CINT'])
})

# This is the regression that protects the augmented approach. A cell that
# references a carrier state IS the individually varying parameter, and the
# subject's last row is the fully informed estimate of it -- not an arbitrary
# evaluation point. Nothing should warn about it, and nobody should later
# "fix" it into being treated like dynamic state dependence.
test_that('carrier state dependence is not a reporting problem', {
  fit <- structure(list(
    standata = list(nlatent = 2L),
    setup = list(matsetup = data.frame(
      parname = c('param*state[9]', 'd22'), row = c(1L, 2L), col = c(1L, 2L),
      matrix = c(7L, 3L)))),
    class = c('ctStanFit', 'ctFit'))

  expect_equal(nrow(ctsem:::.ctFitContextDependentCells(fit)), 1)
  expect_equal(nrow(ctsem:::.ctFitConditionalCells(fit)), 0)
  expect_null(ctsem:::.ctContextNote(ctsem:::.ctFitConditionalCells(fit), 'anywhere'))
})

test_that('a linear model reports nothing at all', {
  fit <- structure(list(
    standata = list(nlatent = 2L),
    setup = list(matsetup = data.frame(parname = c('d11', 'd22'), row = 1:2, col = 1:2,
      matrix = c(3L, 3L)))),
    class = c('ctStanFit', 'ctFit'))
  expect_equal(nrow(ctsem:::.ctFitContextDependentCells(fit)), 0)
})

test_that('detection works from an unfitted model, in both spec syntaxes', {
  nonlinear <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 1,
    manifestNames = 'Y1', latentNames = c('eta1', 'eta2'),
    LAMBDA = matrix(c(1, 0), 1, 2),
    DRIFT = matrix(c('-2*log1p(exp(-2*eta2))', 0, 0, -.00001), 2, 2),
    DIFFUSION = matrix(c('diff', 0, 0, 0), 2, 2)))
  cells <- ctsem:::.ctFitContextDependentCells(nonlinear)
  expect_equal(cells$kind, 'state')
  expect_equal(cells$matrix, 'DRIFT')
  expect_true(ctModelIsNonlinear(nonlinear))

  # TD predictor dependence is the same problem and must be caught the same way.
  tddependent <- suppressMessages(ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
    n.TDpred = 1, manifestNames = 'Y1', latentNames = 'eta1', TDpredNames = 'TD1',
    LAMBDA = matrix(1, 1, 1),
    DRIFT = matrix('-log1p(exp(dr11)) * (1+TD1)', 1, 1),
    DIFFUSION = matrix('diff', 1, 1)))
  expect_equal(ctsem:::.ctFitContextDependentCells(tddependent)$kind, 'tdpred')
  expect_true(ctModelIsNonlinear(tddependent))

  linear <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    LAMBDA = diag(2)))
  expect_false(ctModelIsNonlinear(linear))
})

test_that('the shared note names the matrices and the evaluation point', {
  cells <- data.frame(matrix = c('DRIFT', 'PARS'), row = 1L, col = 1L,
    kind = 'state', stringsAsFactors = FALSE)
  note <- ctsem:::.ctContextNote(cells, 'the T0MEANS state', 'Use state= for another.')
  expect_match(note, 'DRIFT')
  expect_match(note, 'the T0MEANS state')
  expect_match(note, 'Use state= for another')
})
