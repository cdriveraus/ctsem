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
# user wrote, not of that. The parameters that vary belong in the subject
# distribution, which is a separate line, and not in the dynamics.
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
  tex <- paste(as.character(suppressMessages(ctModelLatex(fit,
    equationonly = TRUE, compile = FALSE, open = FALSE, tex = FALSE))),
    collapse = '')

  # Two latent processes, not four: no carrier state reaches the equations.
  expect_false(grepl('eta_{3}', tex, fixed = TRUE))
  # The varying parameters are named once, in the subject distribution.
  expect_true(grepl('vect{\\phi}(i)', tex, fixed = TRUE))
  expect_equal(lengths(regmatches(tex, gregexpr('cint1', tex, fixed = TRUE))), 1L)
})


# The subject parameter distribution and the covariate effects, for a julia
# fit. Not a comparison against stan for the numbers, deliberately: the raw
# population covariance of a weakly identified random effect is where the two
# optimisers most easily land in different places, and a test that needs them
# to agree there is testing the optimisers. The population is set to a known
# value instead, and the equation has to show that value.
test_that("a julia fit's subject distribution is the one it holds", {
  skip_on_cran()
  skip_without_julia()

  generating <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1','Y2'), latentNames = c('L1','L2'), LAMBDA = diag(2),
    DRIFT = matrix(c(-.5,.1,0,-.3), 2, 2, byrow = TRUE),
    DIFFUSION = matrix(c(.8,0,.2,.6), 2, 2, byrow = TRUE),
    MANIFESTVAR = diag(.3,2), MANIFESTMEANS = matrix(0,2,1),
    CINT = matrix(c(.2,-.1),2,1), T0MEANS = matrix(0,2,1), T0VAR = diag(1,2)))
  set.seed(7)
  datalong <- as.data.frame(suppressMessages(ctGenerate(generating,
    n.subjects = 20, Tpoints = 10, burnin = 5, dtmean = 1, logdtsd = 0,
    wide = FALSE)))
  ids <- unique(datalong$id)
  set.seed(9)
  datalong$Z1 <- rnorm(length(ids))[match(datalong$id, ids)]
  datalong$Z2 <- rnorm(length(ids))[match(datalong$id, ids)]

  model <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1','Y2'), latentNames = c('L1','L2'), LAMBDA = diag(2),
    n.TIpred = 2, TIpredNames = c('Z1','Z2'),
    MANIFESTMEANS = matrix(0,2,1), CINT = matrix(c('cint1','cint2'),2,1),
    T0MEANS = matrix(0,2,1), T0VAR = diag(1,2)))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$matrix %in% 'CINT'] <- TRUE

  fit <- suppressWarnings(suppressMessages(ctFit(datalong, model,
    backend = 'julia', cores = 1, optimcontrol = list(estonly = TRUE))))

  # A population the data need not have located: what is under test is whether
  # the equation reports the covariance the fit holds, on the right scale.
  effects <- as.data.frame(ctsem:::.ctBackendSpec(fit)$random_effects)
  raw <- ctsem:::.ctFitRawEstimate(fit)
  raw[effects$parameter[effects$type %in% 'sd']] <- c(-0.5, 0.3)
  raw[effects$parameter[effects$type %in% 'correlation']] <- 0.4
  fit$estimate$raw <- raw
  fit$estimate$rawposterior <- NULL
  fit$transformedpars <- NULL

  # The covariance is the engine's own T0 covariance of the carrier states,
  # which for this representation *is* the population covariance -- so the
  # accessor is checked against the filter rather than against a second copy
  # of the same arithmetic.
  popcov <- ctsem:::.ctBackendRawPopCov(fit)[[1]]$cov
  carrier <- ctCollapse(ctExtract(fit)$pop_T0cov, 1, mean)[3:4, 3:4]
  expect_equal(unname(popcov), unname(carrier))
  expect_true(all(diag(popcov) > 0))

  # A covariate coefficient is a free parameter; `coefficient` is where it sits
  # in the raw vector.
  effects <- as.data.frame(ctsem:::.ctBackendSpec(fit)$ti_effects)
  timat <- ctsem:::.ctBackendRawTipredEffects(fit)
  expect_equal(dim(timat), c(11L, 2L))
  expect_equal(timat[cbind(match(effects$parameter, sort(unique(effects$parameter))),
    effects$predictor)], unname(raw[effects$coefficient]))

  tex <- paste(as.character(suppressMessages(ctModelLatex(fit,
    equationonly = TRUE, compile = FALSE, open = FALSE, tex = FALSE,
    digits = 3))), collapse = '')
  for (cell in unique(as.character(round(c(popcov[1,1], popcov[2,2], popcov[2,1]), 3)))) {
    expect_true(grepl(cell, tex, fixed = TRUE))
  }
  expect_true(grepl('cint1', tex, fixed = TRUE))
  expect_true(grepl('\\text{Z1}', tex, fixed = TRUE))
  expect_true(grepl(as.character(round(timat[1,1], 3)), tex, fixed = TRUE))
})

# `linearise` shows the distribution of the transformed parameters, and this
# backend has no covariance on that scale to show. Refused by name rather than
# accepted and quietly ignored.
test_that("linearise is refused for a julia fit rather than ignored", {
  skip_on_cran()
  m <- suppressMessages(ctModel(type = 'ct', manifestNames = 'Y1', LAMBDA = diag(1)))
  fake <- structure(list(model = m), class = c('ctJuliaFit', 'ctFit'))
  expect_error(ctModelLatex(fake, linearise = TRUE, compile = FALSE, open = FALSE),
    'linearise = TRUE is not available')
  # And the default is the scale this backend does have.
  expect_equal(formals(ctModelLatex)$linearise, quote(inherits(x, 'ctStanFit')))
})

# The union of the covariance, the covariate effects and the initial state is
# ordered by where a parameter first appears, which is not parameter order.
# Both halves have to carry names for the equation to line up, and for a while
# neither did: the raw covariance was dropped and the linearised means were
# rendered against the wrong labels.
test_that("population means and covariance are matched by name, not position", {
  skip_on_cran()
  # The shape the union above produces: parameters ordered by where they first
  # appear -- the covariance first, then whatever only the covariates touch --
  # while the means arrive in parameter order. Position matching sends every
  # mean to the wrong label; both halves carry names so that it cannot.
  pars <- c('b', 'c', 'a')
  popcov <- matrix(0, 3, 3, dimnames = list(pars, pars))
  popcov['b','b'] <- 4; popcov['c','c'] <- 9; popcov['b','c'] <- popcov['c','b'] <- 1
  timat <- matrix(1:6, 3, 2, dimnames = list(pars, c('Z1','Z2')))
  popmeans <- c(a = 10, b = 20, c = 30)
  m <- suppressMessages(ctModel(type = 'ct', manifestNames = 'Y1', LAMBDA = diag(1)))
  m <- c(m, ctsem:::listOfMatrices(m$pars))
  out <- ctsem:::ctModelLatexAugmentT0(popmeans = popmeans, popcov = popcov,
    timat = timat, ctm = m, digits = 3)
  # Character, because the initial state it merges in is the model's own
  # symbolic T0VAR; the point here is which label each number lands on.
  expect_equal(as.numeric(out$popmeans[c('a','b','c')]), c(10, 20, 30))
  expect_equal(as.numeric(out$popcov['b','c']), 1)
})

# The label a population mean is printed under has to be its own.
#
# The parameters are listed in the order they first appear -- the covariance
# first, then whatever only the covariates touch -- and the means arrive in
# parameter order, so the two orders differ as soon as a model has both random
# effects and covariates. Both halves have to carry names. Neither did: the
# linearised means were rendered against the wrong labels, and the raw
# covariance was dropped for zeros.
test_that("each population mean is printed under its own parameter", {
  skip_on_cran()

  generating <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1','Y2'), latentNames = c('L1','L2'), LAMBDA = diag(2),
    DRIFT = matrix(c(-.5,.1,0,-.3), 2, 2, byrow = TRUE),
    DIFFUSION = matrix(c(.8,0,.2,.6), 2, 2, byrow = TRUE),
    MANIFESTVAR = diag(.3,2), MANIFESTMEANS = matrix(0,2,1),
    CINT = matrix(c(.2,-.1),2,1), T0MEANS = matrix(0,2,1), T0VAR = diag(1,2)))
  set.seed(7)
  datalong <- as.data.frame(suppressMessages(ctGenerate(generating,
    n.subjects = 20, Tpoints = 10, burnin = 5, dtmean = 1, logdtsd = 0,
    wide = FALSE)))
  ids <- unique(datalong$id)
  set.seed(9)
  datalong$Z1 <- rnorm(length(ids))[match(datalong$id, ids)]

  model <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1','Y2'), latentNames = c('L1','L2'), LAMBDA = diag(2),
    n.TIpred = 1, TIpredNames = 'Z1',
    MANIFESTMEANS = matrix(0,2,1), CINT = matrix(c('cint1','cint2'),2,1),
    T0MEANS = matrix(0,2,1), T0VAR = diag(1,2)))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$matrix %in% 'CINT'] <- TRUE

  fit <- suppressWarnings(suppressMessages(ctFit(datalong, model,
    backend = 'stan', cores = 1)))

  # The two vectors of the subject distribution line, in the order they are
  # printed: the labels, then the means.
  open <- '\\begin{bmatrix}'
  close <- '\\end{bmatrix}'
  block <- function(tex, from) {
    at <- regexpr(open, substring(tex, from), fixed = TRUE) + from - 1L
    to <- regexpr(close, substring(tex, at), fixed = TRUE) + at - 1L
    cells <- strsplit(substring(tex, at + nchar(open), to - 1L), '\\\\',
      fixed = TRUE)[[1]]
    trimws(cells[nzchar(trimws(cells))])
  }

  for (linearise in c(TRUE, FALSE)) {
    tex <- paste(as.character(suppressMessages(ctModelLatex(fit,
      equationonly = TRUE, compile = FALSE, open = FALSE, tex = FALSE,
      digits = 4, linearise = linearise))), collapse = '')
    labels <- block(tex, 1L)
    means <- block(tex, regexpr('\\mathrm{N} \\left(', tex, fixed = TRUE))
    labels <- sub('_i$', '', sub('\\text{', '', labels, fixed = TRUE))
    labels <- gsub('\\_', '_', sub('}$', '', labels), fixed = TRUE)
    expect_equal(length(labels), length(means))

    ms <- ctsem:::ctMatsetupFreePars(fit$setup$matsetup)
    keep <- as.logical(ms$indvarying + ms$tipred)
    e <- ctExtract(fit)
    truth <- round(ctCollapse(if (linearise) e$popmeans else e$rawpopmeans,
      1, mean), 4)[keep]
    names(truth) <- ms$parname[keep]

    shown <- suppressWarnings(as.numeric(means))
    named <- labels %in% names(truth)
    expect_true(any(named))
    expect_equal(shown[named], as.numeric(truth[labels[named]]))
  }
})
