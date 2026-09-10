# A PARS reference inside T0MEANS or T0VAR.
#
# Declaring a parameter in PARS and naming it in several cells is how ctsem
# shares one parameter between them; `.ctModelIntOverPop()` rewrites each use
# into `PARS[k,1]`. In T0MEANS and T0VAR -- the matrices the engine evaluates
# once, from the parameters -- the julia backend used to refuse every such
# cell, reporting all of them as "does not support a state-dependent
# expression" whether or not any state was involved. That stopped
# `test_behavGenNLcor.R`'s model, whose T0MEANS is four cells sharing two PARS
# parameters and contains no state reference at all.
#
# It is resolved at model-build time now. What is still refused is what
# genuinely cannot be resolved there, and the three shapes below are the
# reasons, each named in its own message.
skip_without_julia()
skip_on_32bit()

.parsref_data <- function(n = 20, seed = 3) {
  set.seed(seed)
  gm <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1','Y2'), latentNames = c('eta1','eta2'), LAMBDA = diag(2),
    DRIFT = diag(-.3, 2), DIFFUSION = diag(.5, 2), MANIFESTVAR = diag(.2, 2),
    MANIFESTMEANS = matrix(0, 2, 1), T0MEANS = matrix(c(1, 1), 2, 1),
    T0VAR = diag(.5, 2), CINT = matrix(0, 2, 1), Tpoints = 8))
  suppressMessages(ctGenerate(gm, n.subjects = n, burnin = 3, dtmean = 1))
}

.parsref_model <- function(T0MEANS = c('baselevel','baselevel'),
  PARS = c('baselevel', 'basecor|2/(1+exp(-param))-1'),
  T0VAR = matrix(c('t0sd', 0, 'basecor', 't0sd'), 2, 2, byrow = TRUE), ...) {
  suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1','Y2'), latentNames = c('eta1','eta2'), LAMBDA = diag(2),
    DRIFT = diag(-.3, 2), DIFFUSION = diag(.5, 2), MANIFESTVAR = diag(.2, 2),
    MANIFESTMEANS = matrix(0, 2, 1), CINT = matrix(0, 2, 1),
    PARS = PARS, T0VAR = T0VAR, T0MEANS = T0MEANS, ...))
}

test_that("a shared PARS parameter in T0MEANS and T0VAR is resolved, not refused", {
  d <- .parsref_data()
  m <- .parsref_model()
  # `intoverpop = TRUE` forces the augmented route, which is what rewrites the
  # named parameters into PARS references. Left at 'auto' this model does not
  # take that route and the references are never created, so the bug does not
  # appear -- which is why a smaller reproduction than this one does not
  # reproduce anything.
  fits <- fit_backends(datalong = d, model = m, cores = 1, verbose = 0,
    intoverpop = TRUE)
  s <- summary(fits$julia, parmatrices = TRUE)
  pm <- s$parmatrices

  t0m <- pm[pm$matrix == 'T0MEANS', ]
  pars <- pm[pm$matrix == 'PARS', ]
  testthat::expect_equal(nrow(t0m), 2L)

  # Both cells name the same parameter, so both carry the same value...
  testthat::expect_equal(t0m$Mean[1], t0m$Mean[2])
  # ...and that value is the referenced PARS cell itself. NOT the cell's own
  # `10 * param` default transform applied on top of it: a reference takes the
  # referenced value, which is what stan does, and composing the outer
  # transform as well would multiply T0MEANS by ten and still converge.
  testthat::expect_equal(t0m$Mean[1], pars$Mean[pars$row == 1], tolerance = 1e-8)

  # Under CTSEM_TEST_STAN, the whole fit against stan. Measured on this
  # design: the two log likelihoods agree exactly (-304.119 on the 25-subject
  # variant used while developing this), and the parameters to about 0.02 --
  # `basecor` is the loosest, being a correlation on a flat ridge at this n.
  expect_backends_agree(fits, tol = 5e-2)
})

test_that("a static cell over two PARS cells is refused, and says so", {
  d <- .parsref_data(n = 12)
  m <- .parsref_model(T0MEANS = c('aa + bb', 'aa'), PARS = c('aa', 'bb'),
    T0VAR = diag(2))
  # Refused rather than mis-resolved, and this is a real constraint rather
  # than an unfinished case: `parameter_transforms.jl` requires each regular
  # transform to read exactly one entry of the parameter vector, and
  # `adjoint_parameters.jl` verifies that when the adjoint workspace is built.
  # A cell over two raw parameters has no single `parnumber` to give.
  testthat::expect_error(
    suppressMessages(ctFit(datalong = d, model = m, fit = FALSE,
      backend = 'julia', intoverpop = TRUE)),
    "combines 2 PARS cells")
})

test_that("a static cell depending on a time-dependent predictor is refused, and says so", {
  d <- .parsref_data(n = 12)
  d <- cbind(d, TD1 = stats::rnorm(nrow(d)))
  m <- .parsref_model(
    T0MEANS = c('baselevel', 'baselevel'),
    PARS = c('baselevel', 'bycov'),
    T0VAR = matrix(c('t0sd', 0, 'bycov * tdpreds[rowi,1]', 't0sd'), 2, 2, byrow = TRUE),
    TDpredNames = 'TD1')
  # This one is a genuine engine limitation rather than a resolution failure:
  # T0VAR is built once from the parameters, and a time-dependent predictor
  # varies by row. Stan evaluates it at the subject's first row; giving the
  # julia t0 block the same access is an engine change, not something that
  # can be composed away here. The message has to say which of the two it is.
  testthat::expect_error(
    suppressMessages(ctFit(datalong = d, model = m, fit = FALSE,
      backend = 'julia', intoverpop = TRUE)),
    "time-dependent predictor data")
})

test_that("a circular t0 reference is still refused", {
  # The thing the blanket refusal was standing in for, and the only thing that
  # actually has no resolution: each cell needs the other's value first.
  # Caught when the model is built, before any backend sees it.
  testthat::expect_error(
    suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
      manifestNames = c('Y1','Y2'), latentNames = c('eta1','eta2'),
      LAMBDA = diag(2), T0MEANS = c('eta2', 'eta1'))),
    "Circular dependency")
})
