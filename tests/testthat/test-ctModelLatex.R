suppressWarnings(suppressPackageStartupMessages(library(ctsem)))
library(testthat)

test_that("ctModelLatex subject distribution excludes calculated pars", {
  m <- suppressWarnings(suppressMessages(ctModel(
    type='ct',
    manifestNames=c('response', 'log_rt'),
    latentNames=c('ddm_drift', 'ddm_boundary_raw', 'ddm_ndt_raw'),
    TDpredNames=c('condition_num', 'stimulus_direction'),
    TDPREDEFFECT=matrix(0, nrow=3, ncol=2),
    LAMBDA=matrix(0, nrow=2, ncol=3),
    manifesttype=c(1, 0),
    MANIFESTMEANS=matrix(c(0, 0), nrow=2),
    MANIFESTVAR=matrix(c(0, 0, 0, 'rt_residual_sd'), nrow=2, byrow=TRUE),
    DRIFT=matrix(0, nrow=3, ncol=3),
    CINT=matrix(0, nrow=3, ncol=1),
    DIFFUSION=matrix(0, nrow=3, ncol=3),
    T0MEANS=matrix(c(
      'drift_t0||TRUE',
      'boundary_raw_t0||TRUE',
      'ndt_raw_t0 | 0.1 + param | TRUE'
    ), nrow=3),
    T0VAR=diag(1e-5, 3),
    PARS=c('gamma_v'),
    silent=TRUE)))
  
  m$pars$param[m$pars$matrix == 'PARS'] <- 'exp(gamma_v)'
  m$pars$indvarying[m$pars$matrix == 'PARS'] <- TRUE
  
  latex <- ctModelLatex(m, tex=FALSE, compile=FALSE, open=FALSE)
  
  expect_match(latex, 'drift\\\\_t0')
  expect_match(latex, 'boundary\\\\_raw\\\\_t0')
  expect_match(latex, 'ndt\\\\_raw\\\\_t0')
  expect_no_match(latex, 'exp\\\\(gamma_v\\\\).*_i')
})

test_that("ctModelLatex handles scalar numeric T0VAR display covariance", {
  m <- suppressMessages(ctModel(type='ct', manifestNames='Y1', LAMBDA=diag(1)))
  expect_error(ctModelLatex(m, tex=FALSE, compile=FALSE, open=FALSE), NA)
})


# A julia fit gets its estimates substituted into the equations, exactly as a
# stan fit does. The check is the two backends against each other: for the same
# model and the same data they converge to the same place, so the documents
# they write are the same document, and any difference is a difference in how
# ctModelLatex reads the two layouts rather than in the fits.
test_that("ctModelLatex substitutes a julia fit's estimates, as stan's", {
  skip_without_julia()

  generating <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1', 'Y2'), latentNames = c('L1', 'L2'),
    LAMBDA = diag(2),
    DRIFT = matrix(c(-.5, .1, 0, -.3), 2, 2, byrow = TRUE),
    DIFFUSION = matrix(c(.8, 0, .2, .6), 2, 2, byrow = TRUE),
    MANIFESTVAR = diag(.3, 2), MANIFESTMEANS = matrix(0, 2, 1),
    CINT = matrix(c(.2, -.1), 2, 1),
    T0MEANS = matrix(0, 2, 1), T0VAR = diag(1, 2)))
  set.seed(7)
  datalong <- suppressMessages(ctGenerate(generating, n.subjects = 30,
    Tpoints = 20, burnin = 5, dtmean = 1, logdtsd = 0, wide = FALSE))

  model <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1', 'Y2'), latentNames = c('L1', 'L2'),
    LAMBDA = diag(2), MANIFESTMEANS = matrix(0, 2, 1),
    CINT = matrix(c('cint1', 'cint2'), 2, 1),
    T0MEANS = matrix(0, 2, 1), T0VAR = diag(1, 2)))
  model$pars$indvarying <- FALSE

  fitjulia <- suppressWarnings(suppressMessages(ctFit(datalong, model,
    backend = 'julia', cores = 1, optimcontrol = list(estonly = TRUE))))
  fitstan <- suppressWarnings(suppressMessages(ctFit(datalong, model,
    backend = 'stan', cores = 1, optimcontrol = list(estonly = TRUE))))

  latex <- function(fit) paste(as.character(suppressMessages(ctModelLatex(fit,
    equationonly = TRUE, compile = FALSE, open = FALSE, tex = FALSE,
    digits = 2))), collapse = '')

  # Estimates, not labels: no free parameter name survives into the equations,
  # and the numbers that replaced them are the ones the fit reports elsewhere.
  expect_false(grepl('cint1', latex(fitjulia), fixed = TRUE))
  drift <- suppressMessages(ctSummaryMatrices(fitjulia)$DRIFT)
  for (cell in as.character(round(diag(drift), 2))) {
    expect_true(grepl(cell, latex(fitjulia), fixed = TRUE))
  }
  # And the same estimates stan writes. Two decimals because the two optimisers
  # stop in slightly different places, not because either is approximate.
  expect_identical(latex(fitjulia), latex(fitstan))
})

# A cell written over the latent state keeps its structure on this backend too:
# the julia parameter table spells the reference `PARS[1,1] * (1 + 0.2 *
# state[2])` where the stan branch reads `dr11 * (1 + 0.2 * eta2)` from the
# base model, and both have to come out as the same equation.
test_that("a state dependent cell keeps its expression for a julia fit", {
  skip_without_julia()

  generating <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
    LAMBDA = diag(2), DRIFT = matrix(c(-.4, .1, 0, -.3), 2, 2),
    CINT = matrix(c(.2, .1), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(.2, 2), DIFFUSION = matrix(c(.5, 0, 0, .4), 2, 2)))
  set.seed(1)
  datalong <- suppressMessages(ctGenerate(generating, n.subjects = 25,
    Tpoints = 12, burnin = 5, dtmean = 1, logdtsd = .1, wide = FALSE))

  model <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
    LAMBDA = diag(2), PARS = c('dr11|-log1p_exp(param)'),
    DRIFT = matrix(c('dr11 * (1 + 0.2 * eta2)', 'd21', 0, 'd22'), 2, 2),
    CINT = matrix(c('c1', 'c2'), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(.2, 2), DIFFUSION = matrix(c('df1', 0, 0, 'df2'), 2, 2)))
  model$pars$indvarying <- FALSE

  fit <- suppressWarnings(suppressMessages(ctFit(datalong, model,
    backend = 'julia', cores = 1, optimcontrol = list(estonly = TRUE))))
  tex <- paste(as.character(suppressMessages(ctModelLatex(fit,
    equationonly = TRUE, compile = FALSE, open = FALSE, tex = FALSE,
    digits = 2))), collapse = '')

  expect_true(grepl('eta_{2}', tex, fixed = TRUE))
  expect_true(grepl('keep their expression', tex, fixed = TRUE))
  # The coordinate the parameter table holds is resolved to the parameter, and
  # then to its estimate, rather than being printed as `PARS[1,1]`.
  expect_false(grepl('PARS[', tex, fixed = TRUE))
  expect_false(grepl('dr11', tex, fixed = TRUE))
})

# An intoverpop fit is written over an augmented state -- one extra latent per
# individually varying parameter -- and the equations are of the system the
# user wrote, not of that. The subject distribution itself this backend does
# not write out, and says so rather than leaving a random effects model looking
# like a fixed effects one.
test_that("an intoverpop julia fit writes the unaugmented system", {
  skip_without_julia()

  set.seed(11)
  n <- 12; tp <- 8
  datalong <- data.frame(id = rep(1:n, each = tp), time = rep(0:(tp - 1), n),
    Y1 = rnorm(n * tp), Y2 = rnorm(n * tp))
  model <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1', 'Y2'), latentNames = c('L1', 'L2'),
    LAMBDA = diag(2), MANIFESTMEANS = matrix(0, 2, 1),
    CINT = matrix(c('cint1', 'cint2'), 2, 1),
    T0MEANS = matrix(0, 2, 1), T0VAR = diag(1, 2)))
  expect_true(sum(model$pars$indvarying) > 0)

  fit <- suppressWarnings(suppressMessages(ctFit(datalong, model,
    backend = 'julia', cores = 1, optimcontrol = list(estonly = TRUE))))
  expect_message(tex <- ctModelLatex(fit, equationonly = TRUE, compile = FALSE,
    open = FALSE, tex = FALSE), 'not written out')
  tex <- paste(as.character(tex), collapse = '')

  # Two latent processes, not four: no carrier state reaches the equations.
  expect_false(grepl('eta_{3}', tex, fixed = TRUE))
  expect_false(grepl('cint1', tex, fixed = TRUE))
})
