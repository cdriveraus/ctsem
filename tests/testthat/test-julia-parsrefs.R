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
# It is resolved at model-build time now, including an expression over
# several PARS cells. What is still refused is what genuinely cannot be
# resolved there -- a row-varying dependence, and a cycle -- each named in its
# own message.
skip_without_julia()
skip_on_32bit()

.parsref_data <- function(n = 20, seed = 3) {
  set.seed(seed)
  gm <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1','Y2'), latentNames = c('eta1','eta2'), LAMBDA = diag(2),
    DRIFT = diag(-.3, 2), DIFFUSION = diag(.5, 2), MANIFESTVAR = diag(.2, 2),
    MANIFESTMEANS = matrix(0, 2, 1), T0MEANS = matrix(c(1, 1), 2, 1),
    T0VAR = diag(.5, 2), CINT = matrix(0, 2, 1), Tpoints = 8))
  suppressMessages(ctGenerate(gm, n = n, burnin = 3, dtmean = 1))
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
  # `intoverpop = 'augmented'` forces the augmented route, which is what rewrites the
  # named parameters into PARS references. Left at 'auto' this model does not
  # take that route and the references are never created, so the bug does not
  # appear -- which is why a smaller reproduction than this one does not
  # reproduce anything.
  fits <- fit_backends(datalong = d, model = m, cores = 1, verbose = 0,
    intoverpop = 'augmented')
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

test_that("a static cell over two PARS cells composes both parameters", {
  # This used to be refused, on the grounds that a regular transform reads one
  # entry of the parameter vector and a cell over two raw parameters has no
  # single `parnumber` to give. The engine's parameter-layer pullback now
  # discovers each transform's support instead of asserting it, so the
  # composition is the same one a single reference gets -- see the comment
  # above `.ctJuliaResolveStaticRefs`.
  #
  # Fit-free, because what it checks is a string: that the rendered transform
  # reads *both* parameters. A composition that dropped one would still fit,
  # still converge, and report a plausible T0MEANS.
  d <- .parsref_data(n = 12)
  m <- .parsref_model(T0MEANS = c('aa + bb', 'aa'), PARS = c('aa', 'bb'),
    T0VAR = diag(2))
  prepped <- suppressMessages(ctFit(datalong = d, model = m, fit = FALSE,
    backend = 'julia', intoverpop = 'augmented'))
  pt <- prepped$parameter_table
  cell <- pt[pt$matrix == 'T0MEANS' & pt$row == 1 & pt$col == 1, ]
  testthat::expect_equal(nrow(cell), 1L)

  reads <- sort(unique(as.integer(regmatches(cell$transform,
    gregexpr('(?<=param\\[)\\d+', cell$transform, perl = TRUE))[[1]])))
  aa <- pt$parnumber[pt$matrix == 'PARS' & pt$row == 1 & pt$col == 1]
  bb <- pt$parnumber[pt$matrix == 'PARS' & pt$row == 2 & pt$col == 1]
  testthat::expect_equal(reads, sort(c(as.integer(aa), as.integer(bb))))
  # And the representative `parnumber` is one of them, which is what the
  # engine's support check requires of an R-rendered transform.
  testthat::expect_true(as.integer(cell$parnumber) %in% reads)
})

test_that("a composed T0VAR fits, and agrees with stan", {
  # The end-to-end half: an identified expression over two parameters, where
  # one of them also appears alone in T0MEANS. Stan writes the expression
  # straight into its model text, so it is the independent implementation to
  # check against rather than a parity mirror.
  d <- .parsref_data(n = 20)
  m <- .parsref_model(T0MEANS = c('aa', 'aa'), PARS = c('aa', 'bb'),
    T0VAR = matrix(c('log1p_exp(aa + bb)', 0, 0, 't0sd'), 2, 2, byrow = TRUE))
  fits <- fit_backends(datalong = d, model = m, cores = 1, verbose = 0,
    intoverpop = 'augmented')
  testthat::expect_s3_class(fits$julia, 'ctJuliaFit')
  expect_backends_agree(fits, tol = 5e-2)
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
      backend = 'julia', intoverpop = 'augmented')),
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

test_that("a per-row matrix's expression that reads no state is evaluated once a subject", {
  # DRIFT, CINT, DIFFUSION, the measurement matrices and the rest are the ones
  # the engine can also evaluate per row, and an expression cell in them used
  # to be, whatever it read. One reading only parameters, other cells and the
  # TI predictors -- the drift-scaled diffusion and intercept reading their
  # DRIFT diagonal, a mean switched by a TI predictor -- is composed into an
  # ordinary transform instead (`.ctJuliaParameterTable()`). The oracle is the
  # same model with each of those cells also reading a state at zero weight,
  # which keeps it per row: the objective and gradient must not move.
  set.seed(4)
  d <- as.data.frame(.parsref_data(n = 10))
  d$z <- as.numeric(d$id %% 2 == 0)
  build <- function(per_row) {
    m <- suppressWarnings(suppressMessages(ctModel(type = 'ct', n.latent = 2,
      n.manifest = 2, manifestNames = c('Y1','Y2'), latentNames = c('eta1','eta2'),
      LAMBDA = diag(2), MANIFESTMEANS = matrix(c('mm', 0), 2),
      CINT = matrix(c('cint_a', 'cint_b'), 2), PARS = matrix(c('p1', 'p2'), 1, 2),
      TIpredNames = 'z', driftScaled = c('DIFFUSION', 'CINT'))))
    m$pars$indvarying <- FALSE
    mm <- m$pars$matrix == 'MANIFESTMEANS' & m$pars$row == 1
    m$pars$param[mm] <- '(1 - z) * p1 + z * p2'
    expression <- !is.na(m$pars$param) & grepl('[', m$pars$param, fixed = TRUE)
    if (per_row) {
      cells <- expression & m$pars$matrix %in% c('DIFFUSION', 'CINT') | mm
      m$pars$param[cells] <- paste0(m$pars$param[cells], ' + 0 * state[1]')
    }
    m
  }
  spec <- lapply(c(once = FALSE, per_row = TRUE), function(per_row)
    suppressWarnings(suppressMessages(ctFit(datalong = d, model = build(per_row),
      fit = FALSE, backend = 'julia'))))
  rowwise <- function(pt) {
    blank <- function(x) ifelse(is.na(x), '', as.character(x))
    nzchar(paste0(blank(pt$predicttransform), blank(pt$updatetransform),
      blank(pt$tdtransform)))
  }
  once <- spec$once$parameter_table
  per_row <- spec$per_row$parameter_table
  cell <- function(pt, mat, r, c) which(pt$matrix == mat & pt$row == r & pt$col == c)
  # No cell of the first is per row; the diffusion diagonals, the intercepts
  # and the switched mean of the second are, and the first reads the
  # predictor from the parameter vector's tail.
  testthat::expect_false(any(rowwise(once)))
  for (rc in list(c('DIFFUSION', 1, 1), c('DIFFUSION', 2, 2), c('CINT', 1, 1),
      c('MANIFESTMEANS', 1, 1))) {
    testthat::expect_true(rowwise(per_row)[cell(per_row, rc[1], rc[2], rc[3])],
      info = paste(rc, collapse = ' '))
    testthat::expect_false(is.na(once$transform[cell(once, rc[1], rc[2], rc[3])]),
      info = paste(rc, collapse = ' '))
  }
  testthat::expect_match(once$transform[cell(once, 'MANIFESTMEANS', 1, 1)], 'tipreds[1]',
    fixed = TRUE)
  # That read is a parameter position by the time the engine builds the
  # transform, which it interprets: no compile to announce for it.
  testthat::expect_false(ctsem:::.ctJuliaAnnounceCompilation(
    ctsem:::.ctJuliaModule(spec$once$project), once))
  npar <- ctsem:::.ctBackendNpar(spec$once)
  testthat::expect_equal(npar, ctsem:::.ctBackendNpar(spec$per_row))
  for (k in 1:2) {
    x <- stats::rnorm(npar, 0, 0.3)
    a <- ctJuliaEvaluate(spec$once, x, cores = 1)
    b <- ctJuliaEvaluate(spec$per_row, x, cores = 1)
    testthat::expect_equal(a$value, b$value, tolerance = 1e-10)
    testthat::expect_equal(as.numeric(a$gradient), as.numeric(b$gradient),
      tolerance = 1e-8)
  }
})
