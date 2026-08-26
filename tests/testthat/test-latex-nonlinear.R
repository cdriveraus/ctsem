library(ctsem)
library(testthat)

context('latexnonlinear')

# Substituting the fitted number into every cell turns a nonlinear model's
# equation into a linear model's equation. These pin the alternative: keep the
# structure, substitute only the labels.

test_that('parameter labels become estimates and references become math', {
  render <- ctsem:::.ctLatexRenderExpression
  estimates <- c(dr11 = -0.3123, a21 = 1.5, a1 = 99)
  latent <- c('eta1', 'eta2')
  tdpred <- 'TD1'

  expect_equal(render('-log1p(exp(dr11)) * eta2', estimates, latent, tdpred, 3),
    '-\\mathrm{log1p}(\\mathrm{exp}(-0.312)) \\cdot \\eta_{2}')
  expect_equal(render('-2*log1p(exp(-2*state[2]))', estimates, latent, tdpred, 3),
    '-2 \\cdot \\mathrm{log1p}(\\mathrm{exp}(-2 \\cdot \\eta_{2}))')
  expect_equal(render('dr11 * (1+TD1)', estimates, latent, tdpred, 3),
    '-0.312 \\cdot (1+\\text{TD1})')
  expect_equal(render('a21 * tdpreds[rowi, 1]', estimates, latent, tdpred, 3),
    '1.5 \\cdot \\text{TD}_{1}')
})

test_that('a longer parameter name is substituted before a shorter prefix of it', {
  # `a1` must not eat the `a1` inside `a21`... nor `a2` inside `a21`.
  expect_equal(ctsem:::.ctLatexRenderExpression('a21', c(a21 = 1.5, a1 = 99),
    character(), character(), 3), '1.5')
})

test_that('a macro followed by a parenthesis is not read as a function name', {
  # `\\cdot (` was previously rewritten to `\\mathrm{cdot} (`.
  out <- ctsem:::.ctLatexRenderExpression('dr11*(1+eta1)', c(dr11 = 2),
    c('eta1'), character(), 3)
  expect_false(grepl('mathrm{cdot}', out, fixed = TRUE))
  expect_true(grepl('\\cdot', out, fixed = TRUE))
})

test_that('the expression is reassembled from either spec syntax', {
  piped <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 1,
    n.TDpred = 1, manifestNames = 'Y1', latentNames = c('eta1', 'eta2'),
    TDpredNames = 'TD1', LAMBDA = matrix(c(1, 0), 1, 2),
    DRIFT = matrix(c('dr11 | -log1p(exp(param)) * (1+TD1)', 0, 0, -.1), 2, 2),
    DIFFUSION = matrix(c('diff', 0, 0, 0), 2, 2)))
  expect_equal(ctsem:::.ctLatexCellExpression(piped$pars, 'DRIFT', 1, 1),
    '-log1p(exp(dr11))*(1+TD1)')

  bare <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 1,
    manifestNames = 'Y1', latentNames = c('eta1', 'eta2'),
    LAMBDA = matrix(c(1, 0), 1, 2),
    DRIFT = matrix(c('-2*log1p(exp(-2*eta2))', 0, 0, -.00001), 2, 2),
    DIFFUSION = matrix(c('diff', 0, 0, 0), 2, 2)))
  expect_equal(ctsem:::.ctLatexCellExpression(bare$pars, 'DRIFT', 1, 1),
    '-2*log1p(exp(-2*eta2))')
})

test_that('a linear fit renders exactly as before, with no caveat', {
  skip_if_not(exists('ctstantestfit'))
  tex <- paste(ctModelLatex(ctstantestfit, equationonly = TRUE, compile = FALSE,
    open = FALSE, tex = FALSE), collapse = '')
  expect_false(grepl('keep their expression', tex, fixed = TRUE))
})

test_that('a context-dependent cell survives into the equation', {
  skip_if_not(exists('ctstantestfit'))
  fit <- ctstantestfit
  pars <- fit$ctstanmodelbase$pars
  cell <- which(pars$matrix == 'DRIFT' & pars$row == 1 & pars$col == 1)
  pars$param[cell] <- '-log1p(exp(drift11)) * eta2'
  pars$transform[cell] <- NA
  fit$ctstanmodelbase$pars <- pars

  tex <- paste(ctModelLatex(fit, equationonly = TRUE, compile = FALSE,
    open = FALSE, tex = FALSE), collapse = '')
  # The point of the whole exercise: the equation still says it depends on eta2.
  expect_true(grepl('eta_{2}', tex, fixed = TRUE))
  expect_true(grepl('keep their expression', tex, fixed = TRUE))
})
