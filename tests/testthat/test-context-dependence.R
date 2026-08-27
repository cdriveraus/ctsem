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

test_that('cells mixing a carrier and a dynamic reference are detected', {
  mixed <- structure(list(
    standata = list(nlatent = 2L),
    setup = list(matsetup = data.frame(parname = c('state[9]*state[1]', 'd22'),
      row = c(1L, 2L), col = c(1L, 2L), matrix = c(3L, 3L)))),
    class = c('ctStanFit', 'ctFit'))
  expect_equal(nrow(ctsem:::.ctContextConflatedCells(mixed)), 1)

  carrier <- structure(list(
    standata = list(nlatent = 2L),
    setup = list(matsetup = data.frame(parname = 'param*state[9]', row = 1L, col = 1L,
      matrix = 7L))),
    class = c('ctStanFit', 'ctFit'))
  expect_equal(nrow(ctsem:::.ctContextConflatedCells(carrier)), 0)
})

test_that('the message fires for dynamic and TD dependence and not otherwise', {
  mk <- function(parname, matrix) structure(list(
    standata = list(nlatent = 2L),
    setup = list(matsetup = data.frame(parname = parname, row = 1L, col = 1L,
      matrix = matrix))),
    class = c('ctStanFit', 'ctFit'))

  expect_message(ctsem:::.ctContextMessage(mk('log1p(exp(state[1]))', 3L),
    ctsem:::.ctContextPopLabel), 'latent state')
  expect_message(ctsem:::.ctContextMessage(mk('b*tdpreds[rowi, 1]', 3L),
    ctsem:::.ctContextPopLabel), 'time dependent predictor')
  expect_silent(ctsem:::.ctContextMessage(mk('param*state[9]', 7L),
    ctsem:::.ctContextPopLabel))
  expect_silent(ctsem:::.ctContextMessage(mk('d11', 3L), ctsem:::.ctContextPopLabel))
})

test_that('the remedy names the julia backend only for a julia fit', {
  cells <- mk <- structure(list(
    standata = list(nlatent = 2L),
    setup = list(matsetup = data.frame(parname = 'log1p(exp(state[1]))', row = 1L,
      col = 1L, matrix = 3L))),
    class = c('ctStanFit', 'ctFit'))
  expect_match(ctsem:::.ctContextRemedy(mk), "backend='julia'")
  class(mk) <- c('ctJuliaFit', 'ctFit')
  expect_match(ctsem:::.ctContextRemedy(mk), 'state=')
})

# The bundled fit is linear, but its individually varying CINT parameters are
# carried as latent states -- so its DRIFT/CINT expressions do contain
# `state[k]`. A detector that looked only for that string would warn about
# every ordinary random-effects model in the package.
test_that('a linear random-effects fit is not reported as nonlinear', {
  skip_if_not(exists('ctstantestfit'))
  cells <- ctsem:::.ctFitContextDependentCells(ctstantestfit)
  expect_true(all(cells$kind == 'carrier'))
  expect_false(ctModelIsNonlinear(ctstantestfit))
  expect_equal(nrow(ctsem:::.ctFitConditionalCells(ctstantestfit)), 0)
  expect_null(summary(ctstantestfit, priorcheck = FALSE)$parmatNote)
})

# Previously a hard stop(), which meant ctDiscretePars(standardise=TRUE) aborted
# on ctsem's own bundled example fit -- a linear model -- whenever a posterior
# draw happened to be non-stationary. For a model whose DRIFT depends on the
# state it aborts far more often, and blames the model rather than the point the
# model was linearised at.
test_that('a non-stationary standardisation returns NaN rather than aborting', {
  pars <- list(
    DRIFT = array(c(.5, 0, 0, .5), dim = c(1, 2, 2)),
    DIFFUSIONcov = array(diag(2), dim = c(1, 2, 2)),
    asymDIFFUSIONcov = array(c(-1, 0, 0, -1), dim = c(1, 2, 2)))
  expect_message(
    out <- ctsem:::ctDiscreteParsDrift(pars, times = 1, observational = FALSE,
      standardise = TRUE, quiet = FALSE),
    'no stationary variance')
  expect_true(all(is.nan(out)))
})

test_that('ctDiscretePars refuses state= for a stan fit, with a reason', {
  skip_if_not(exists('ctstantestfit'))
  expect_error(ctDiscretePars(ctstantestfit, times = 1, state = 'mean'),
    "backend='julia'")
})

# Phase portrait, on the linear side where no engine is needed and the answers
# are all in closed form.
test_that('the phase portrait field and fixed point agree with each other', {
  skip_if_not(exists('ctstantestfit'))
  portrait <- ctPhasePortrait(ctstantestfit, gridsize = 7, plot = FALSE)
  expect_true(all(is.finite(portrait$field$dx)))
  expect_true(all(is.finite(portrait$field$dy)))

  # The marked fixed point must be where the drawn field actually vanishes --
  # not merely near it. Reading asymCINT from the summary does not satisfy
  # this, because the summary collapses each matrix over draws separately.
  field <- ctsem:::.ctFieldFunction(ctstantestfit)
  point <- as.numeric(unlist(portrait$fixedpoint))
  expect_lt(max(abs(field(point))), 1e-8)

  # A linear model's nullclines are straight, so each is a single piece.
  expect_equal(length(unique(portrait$nullclines$piece)), 2)
})

test_that('the phase portrait refuses a nonlinear stan fit', {
  skip_if_not(exists('ctstantestfit'))
  fit <- ctstantestfit
  setup <- fit$setup$matsetup
  cell <- setup$matrix == 3 & setup$row == 1 & setup$col == 1
  fit$setup$matsetup$parname[cell] <- 'log1p(exp(state[1]))'
  expect_true(ctModelIsNonlinear(fit))
  expect_error(ctPhasePortrait(fit, gridsize = 3), "backend='julia'")
})

test_that('the phase portrait validates its latents argument', {
  skip_if_not(exists('ctstantestfit'))
  expect_error(ctPhasePortrait(ctstantestfit, latents = 1), 'exactly two')
  expect_error(ctPhasePortrait(ctstantestfit, latents = c('nope', 'eta2')),
    'exactly two')
})


# Every interpretation of "a one unit change in process c" is dtDRIFT %*% C for
# a different companion matrix C. These are pure linear algebra and need no
# fit, which is the point: the whole family can be pinned in a couple of
# seconds rather than by refitting models.
test_that('each companion matrix is what it claims to be', {
  drift <- matrix(c(-.4, .15, .05, -.3), 2, 2)
  diffusion <- matrix(c(1, 1.6, 1.6, 16), 2, 2)   # innovation sds 1 and 4
  stationary <- matrix(solve(kronecker(diag(2), drift) + kronecker(drift, diag(2)),
    -as.vector(diffusion)), 2, 2)
  companion <- function(type) ctsem:::.ctCompanionMatrix(type, diffusion,
    stationary, 2)

  # experimental: nothing else moves.
  expect_equal(companion('experimental'), diag(2))
  # observational: the conditional expectation under the STATE covariance.
  expect_equal(companion('observational'),
    stationary %*% diag(1 / diag(stationary)))
  # shock: the conditional expectation under the INNOVATION covariance.
  expect_equal(companion('shock'), diffusion %*% diag(1 / diag(diffusion)))
  # These are different questions and must give different answers whenever the
  # state and innovation correlations differ.
  expect_false(isTRUE(all.equal(companion('observational'), companion('shock'))))

  # Every one has a unit diagonal -- "process c moves by one unit" -- and is
  # asymmetric, because E[x_2|x_1=1] and E[x_1|x_2=1] are different numbers
  # unless the variances match. A correlation matrix is symmetric and so can
  # never be any of these.
  for (type in c('observational', 'shock')) {
    expect_equal(diag(companion(type)), c(1, 1))
    expect_false(isTRUE(all.equal(companion(type), t(companion(type)))))
  }

  # orthogonal: the Cholesky factor, so the implied shock covariance recovers
  # the diffusion correlation.
  orthogonal <- companion('orthogonal')
  scale <- diag(diag(chol(diffusion))^2)
  expect_equal(cov2cor(orthogonal %*% scale %*% t(orthogonal)),
    cov2cor(diffusion), tolerance = 1e-8)

  expect_error(ctsem:::.ctCompanionMatrix('nonsense', diffusion, stationary, 2))
})

test_that('the historical logical argument still selects the right two', {
  expect_equal(ctsem:::.ctCompanionType(FALSE), 'experimental')
  expect_equal(ctsem:::.ctCompanionType(TRUE), 'observational')
  expect_equal(ctsem:::.ctCompanionType('shock'), 'shock')
})

test_that('ctDiscreteParsDrift applies the companion matrix it was asked for', {
  drift <- matrix(c(-.4, .15, .05, -.3), 2, 2)
  diffusion <- matrix(c(1, 1.6, 1.6, 16), 2, 2)
  stationary <- matrix(solve(kronecker(diag(2), drift) + kronecker(drift, diag(2)),
    -as.vector(diffusion)), 2, 2)
  pars <- list(DRIFT = array(drift, dim = c(1, 2, 2)),
    DIFFUSIONcov = array(diffusion, dim = c(1, 2, 2)),
    asymDIFFUSIONcov = array(stationary, dim = c(1, 2, 2)))
  scales <- diag(sqrt(diag(stationary)))
  transition <- as.matrix(Matrix::expm(drift * 1.5))

  for (type in ctsem:::.ctCompanionTypes) {
    got <- ctsem:::ctDiscreteParsDrift(pars, times = 1.5, observational = type,
      standardise = FALSE, quiet = TRUE)
    expect_equal(got[1, 1, 1, , ],
      transition %*% ctsem:::.ctCompanionMatrix(type, diffusion, stationary, 2),
      ignore_attr = TRUE)
  }

  # Standardising composes on top: S^-1 (dtA C) S. For 'observational' that is
  # the standardised simple regression, S^-1 dtA S R.
  std <- ctsem:::ctDiscreteParsDrift(pars, times = 1.5, observational = TRUE,
    standardise = TRUE, quiet = TRUE)
  expect_equal(std[1, 1, 1, , ],
    solve(scales) %*% transition %*% scales %*% cov2cor(stationary),
    ignore_attr = TRUE)
})

test_that('a process with no diffusion gets no companions rather than NaN', {
  pars <- list(DRIFT = array(matrix(c(-.4, .1, 0, -.3), 2, 2), dim = c(1, 2, 2)),
    DIFFUSIONcov = array(matrix(c(1, 0, 0, 0), 2, 2), dim = c(1, 2, 2)),
    asymDIFFUSIONcov = array(matrix(c(1, .5, .5, 4), 2, 2), dim = c(1, 2, 2)))
  out <- ctsem:::ctDiscreteParsDrift(pars, times = 1, observational = 'shock',
    standardise = FALSE, quiet = TRUE)
  expect_true(all(is.finite(out)))
})

# observational + standardise is the model implied cross-correlation function,
# Cor(x_r(t+u), x_c(t)) -- the counterpart to what ctACF computes from data.
# Verified against a simulated series separately; pinned here algebraically so
# it runs in milliseconds.
test_that('observational and standardise together give the cross-correlation', {
  drift <- matrix(c(-.4, .15, .05, -.3), 2, 2)
  diffusion <- matrix(c(1, 1.6, 1.6, 16), 2, 2)
  stationary <- matrix(solve(kronecker(diag(2), drift) + kronecker(drift, diag(2)),
    -as.vector(diffusion)), 2, 2)
  pars <- list(DRIFT = array(drift, dim = c(1, 2, 2)),
    DIFFUSIONcov = array(diffusion, dim = c(1, 2, 2)),
    asymDIFFUSIONcov = array(stationary, dim = c(1, 2, 2)))

  for (lag in c(0, 1.5, 4)) {
    got <- ctsem:::ctDiscreteParsDrift(pars, times = lag, observational = TRUE,
      standardise = TRUE, quiet = TRUE)[1, 1, 1, , ]
    # Cor(x_r(t+u), x_c(t)) = [dtA %*% Sigma]_rc / (sd_r sd_c)
    sdv <- sqrt(diag(stationary))
    expected <- (as.matrix(Matrix::expm(drift * lag)) %*% stationary) /
      outer(sdv, sdv)
    expect_equal(got, expected, ignore_attr = TRUE)
  }

  # At lag zero that is just the latent correlation matrix.
  at0 <- ctsem:::ctDiscreteParsDrift(pars, times = 0, observational = TRUE,
    standardise = TRUE, quiet = TRUE)[1, 1, 1, , ]
  expect_equal(at0, cov2cor(stationary), ignore_attr = TRUE)
})
