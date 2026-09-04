suppressWarnings(suppressPackageStartupMessages(library(ctsem)))
library(testthat)

# The load-bearing claims of ctNetwork() are three closed forms, each checkable
# by hand rather than against another implementation of itself:
#
#   temporal        expm(DRIFT * dt), which for a lower triangular DRIFT has an
#                   elementary closed form.
#   innovation      int_0^dt expm(DRIFT s) Q expm(DRIFT s)' ds, which for
#                   DRIFT = -I is Q (1 - exp(-2 dt)) / 2.
#   contemporaneous partial correlations of that, which for two variables are
#                   just the correlation.
#
# The structural claims -- edge direction, thresholding, what is in the returned
# object -- are checked separately, because getting the numbers right and the
# from/to convention backwards would be a silent disaster in a figure.

.netmodel <- function()
  suppressMessages(suppressWarnings(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('y1', 'y2'), latentNames = c('a', 'b'),
    LAMBDA = diag(2),
    DRIFT = matrix(c(-1, 0, .5, -2), 2, 2, byrow = TRUE),
    DIFFUSION = matrix(c(1, 0, .3, 1), 2, 2, byrow = TRUE))))


test_that('the innovation covariance matches its closed form', {
  q <- matrix(c(1, .5, .5, 1), 2, 2)
  # DRIFT = -I: expm(-I s) = exp(-s) I, so the integral is Q int exp(-2s) ds.
  for(dt in c(.3, .7, 4)){
    expect_equal(ctsem:::.ctNetworkInnovation(-diag(2), q, dt, TRUE),
      q * (1 - exp(-2 * dt)) / 2)
  }

  # For a general stable DRIFT, against the stationary identity
  # Q(dt) = asym - dtDRIFT asym dtDRIFT', which is a wholly different route.
  drift <- matrix(c(-1, .3, -.2, -2), 2, 2)
  asym <- ctsem:::.ctNetworkAsym(drift, q, TRUE)
  expect_lt(max(abs(drift %*% asym + asym %*% t(drift) + q)), 1e-10)
  dtdrift <- as.matrix(expm::expm(drift * 1.3))
  expect_equal(ctsem:::.ctNetworkInnovation(drift, q, 1.3, TRUE),
    asym - dtdrift %*% asym %*% t(dtdrift), tolerance = 1e-10)

  # DRIFT = -I again: the stationary covariance is Q / 2.
  expect_equal(ctsem:::.ctNetworkAsym(-diag(2), q, TRUE), q / 2)

  # A non-stationary system has no stationary covariance to report.
  expect_null(ctsem:::.ctNetworkAsym(matrix(c(.5, 0, 0, -1), 2, 2), q, TRUE))
})


test_that('partial correlations are the Gaussian graphical model', {
  # Two variables: the partial correlation is the correlation.
  expect_equal(ctsem:::.ctNetworkPcor(matrix(c(1, .5, .5, 1), 2, 2))[1, 2], .5)
  expect_equal(diag(ctsem:::.ctNetworkPcor(matrix(c(4, 1, 1, 9), 2, 2))), c(0, 0))

  # Three variables, a chain: 1 and 3 are correlated only through 2, so their
  # partial correlation is zero and the other two are the conditionals.
  prec <- matrix(c(1, -.4, 0, -.4, 1, -.5, 0, -.5, 1), 3, 3)
  p <- ctsem:::.ctNetworkPcor(solve(prec))
  expect_lt(abs(p[1, 3]), 1e-10)
  expect_equal(p[1, 2], .4)
  expect_equal(p[2, 3], .5)

  # A noiseless process makes the covariance singular; it gets no edges rather
  # than taking the whole matrix down.
  p <- ctsem:::.ctNetworkPcor(matrix(c(1, .5, 0, .5, 1, 0, 0, 0, 0), 3, 3))
  expect_equal(p[3, ], c(0, 0, 0))
  expect_equal(p[1, 2], .5)
})


test_that('the temporal network from a model specification is expm(DRIFT*dt)', {
  m <- .netmodel()
  for(dt in c(.5, 1, 3)){
    net <- ctNetwork(m, dt = dt, standardise = FALSE, quiet = TRUE)
    # DRIFT = [[-1,0],[.5,-2]] is lower triangular, so expm(DRIFT t) has
    # exp(-t) and exp(-2t) on the diagonal and .5 (exp(-t) - exp(-2t)) / (2-1)
    # below it.
    expect_equal(unname(diag(net$temporal)), c(exp(-dt), exp(-2 * dt)))
    expect_equal(unname(net$temporal[2, 1]), .5 * (exp(-dt) - exp(-2 * dt)))
    expect_equal(unname(net$temporal[1, 2]), 0)
    expect_equal(dimnames(net$temporal), list(c('a', 'b'), c('a', 'b')))
  }

  # standardise=TRUE puts an edge in sd units of the two processes: the
  # unstandardised weight times sd(from) / sd(to).
  net <- ctNetwork(m, dt = 1, standardise = TRUE, quiet = TRUE)
  sdv <- sqrt(diag(net$asymDIFFUSIONcov))
  expect_equal(unname(net$temporal[2, 1]),
    unname(.5 * (exp(-1) - exp(-2)) * sdv[1] / sdv[2]))
  expect_equal(unname(diag(net$temporal)), c(exp(-1), exp(-2)))
})


test_that('the contemporaneous network is a dt-free property of DIFFUSION when the drift is isotropic', {
  # expm(-k I s) is a scalar times the identity, so Q(dt) is DIFFUSIONcov times
  # a positive scalar and its partial correlations do not depend on dt at all.
  m <- suppressMessages(suppressWarnings(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('y1', 'y2'), latentNames = c('a', 'b'), LAMBDA = diag(2),
    DRIFT = matrix(c(-1, 0, 0, -1), 2, 2),
    DIFFUSION = matrix(c(1, 0, .4, 1), 2, 2, byrow = TRUE))))
  short <- ctNetwork(m, dt = .05, standardise = FALSE, quiet = TRUE)
  long <- ctNetwork(m, dt = 20, standardise = FALSE, quiet = TRUE)
  expect_equal(short$contemporaneous, long$contemporaneous)
  # And with two processes it is the correlation implied by the DIFFUSION spec.
  expect_equal(unname(short$contemporaneous[1, 2]),
    unname(cov2cor(short$DIFFUSIONcov)[1, 2]))
  # Symmetric, zero diagonal, in [-1,1] -- a partial correlation matrix.
  expect_equal(short$contemporaneous, t(short$contemporaneous))
  expect_equal(unname(diag(short$contemporaneous)), c(0, 0))
  expect_true(all(abs(short$contemporaneous) <= 1))
})


test_that('the returned object has the documented structure and edge directions', {
  m <- .netmodel()
  net <- ctNetwork(m, dt = 1, standardise = FALSE, quiet = TRUE,
    networks = c('temporal', 'contemporaneous', 'asymptotic', 'measurement'))

  expect_s3_class(net, 'ctNetwork')
  expect_true(all(c('temporal', 'contemporaneous', 'asymptotic', 'measurement',
    'DRIFT', 'DIFFUSIONcov', 'asymDIFFUSIONcov', 'innovation', 'edges',
    'nodes') %in% names(net)))
  expect_equal(dim(net$temporal), c(2, 2))
  expect_equal(dim(net$measurement), c(2, 2))
  expect_equal(dimnames(net$measurement), list(c('y1', 'y2'), c('a', 'b')))
  expect_equal(net$nodes$name, c('a', 'b', 'y1', 'y2'))
  expect_equal(net$nodes$type, c('latent', 'latent', 'manifest', 'manifest'))
  expect_equal(attr(net, 'dt'), 1)

  expect_named(net$edges,
    c('network', 'from', 'to', 'weight', 'directed', 'self'))
  # mat[i,j] is the effect of column j on row i, so the only cross-lagged edge,
  # temporal[2,1], must run from a to b and not the other way.
  cross <- net$edges[net$edges$network == 'temporal' & !net$edges$self, ]
  expect_equal(nrow(cross), 1)
  expect_equal(cross$from, 'a')
  expect_equal(cross$to, 'b')
  expect_equal(cross$weight, .5 * (exp(-1) - exp(-2)))
  expect_true(all(net$edges$directed[net$edges$network == 'temporal']))

  # An undirected network contributes one edge per pair, not two.
  expect_equal(sum(net$edges$network == 'contemporaneous'), 1)
  expect_false(any(net$edges$directed[net$edges$network == 'contemporaneous']))
  expect_equal(sum(net$edges$network == 'measurement'), 2)

  # Thresholding drops edges from the list but never touches the matrices.
  thin <- ctNetwork(m, dt = 1, standardise = FALSE, quiet = TRUE, threshold = .2)
  expect_equal(thin$temporal, net$temporal)
  expect_false(any(abs(thin$edges$weight) <= .2))
  expect_lt(nrow(thin$edges), nrow(net$edges))

  # networks= restricts the edge list only.
  only <- ctNetwork(m, dt = 1, standardise = FALSE, quiet = TRUE, networks = 'temporal')
  expect_equal(unique(only$edges$network), 'temporal')
  expect_false(is.null(only$contemporaneous))

  expect_error(ctNetwork(m, dt = c(1, 2)), 'single positive number')
  expect_error(ctNetwork('not a model'), 'needs a ctModel')
})


test_that('a discrete time model uses powers of DRIFT and whole steps', {
  m <- suppressMessages(suppressWarnings(ctModel(type = 'dt', n.latent = 2, n.manifest = 2,
    manifestNames = c('y1', 'y2'), latentNames = c('a', 'b'), LAMBDA = diag(2),
    DRIFT = matrix(c(.5, 0, .2, .3), 2, 2, byrow = TRUE),
    DIFFUSION = matrix(c(1, 0, 0, 1), 2, 2))))
  drift <- matrix(c(.5, 0, .2, .3), 2, 2, byrow = TRUE)
  net <- ctNetwork(m, dt = 3, standardise = FALSE, quiet = TRUE)
  expect_equal(unname(net$temporal), drift %*% drift %*% drift)
  expect_error(ctNetwork(m, dt = 1.5, quiet = TRUE), 'whole number of steps')
})


test_that('a fit gives the networks its own reported matrices imply', {
  utils::data('ctstantestfit', package = 'ctsem', envir = environment())
  net <- suppressMessages(ctNetwork(ctstantestfit, dt = 1.4,
    standardise = FALSE, networks = c('temporal', 'contemporaneous', 'measurement')))

  expect_s3_class(net, 'ctNetwork')
  expect_equal(attr(net, 'source'), 'fit')
  mats <- suppressMessages(ctSummaryMatrices(ctstantestfit))
  nl <- nrow(net$temporal)
  expect_equal(net$DRIFT, as.matrix(mats$DRIFT)[1:nl, 1:nl],
    ignore_attr = TRUE)
  # The temporal network is expm of the fit's own DRIFT, computed here without
  # any of ctNetwork's machinery.
  expect_equal(unname(net$temporal),
    unname(as.matrix(expm::expm(net$DRIFT * 1.4))), tolerance = 1e-8)
  expect_equal(unname(net$contemporaneous),
    unname(ctsem:::.ctNetworkPcor(ctsem:::.ctNetworkInnovation(net$DRIFT,
      net$DIFFUSIONcov, 1.4, TRUE))))
  expect_equal(rownames(net$temporal), ctstantestfit$ctstanmodel$latentNames[1:nl])
  expect_true(nrow(net$edges) > 0)
  expect_setequal(unique(net$edges$network),
    c('temporal', 'contemporaneous', 'measurement'))
})


test_that('plotting returns a ggplot, one panel per network and interval', {
  m <- .netmodel()
  g <- suppressMessages(ctNetworkPlot(m, dt = c(.5, 2), networks = 'temporal'))
  expect_s3_class(g, 'ggplot')
  expect_equal(nlevels(g$layers[[1]]$data$panel), 2)

  net <- ctNetwork(m, dt = 1, quiet = TRUE)
  expect_s3_class(suppressMessages(ctNetworkPlot(net)), 'ggplot')
  expect_error(suppressMessages(ctNetworkPlot(m, dt = 1, threshold = 10)),
    'nothing to draw')

  # print() says what it is without erroring on either source.
  expect_output(print(net), 'ctNetwork from a model')
})


# State dependence -------------------------------------------------------------
#
# For a model whose DRIFT depends on a latent process there is no such thing as
# *the* network: what a fit reports is one linearisation, and the graph drawn
# from it is that linearisation's graph. The failure this guards against is not
# a wrong number, it is a right number with no label -- a figure that leaves the
# session claiming to be the system's network.
#
# The detector itself is tested against real specifications in
# test-context-dependence.R. These tests are about the wiring: that ctNetwork
# asks, says so once, names the point in the package's own words, puts it in the
# figure, and stays completely silent on a linear model.

# A fit whose DRIFT[1,1] is written by an expression referencing state[2].
# Mutating the bundled fit's parameter table rather than optimising a nonlinear
# model keeps this in the fast tier: every reporting path reads the table, so
# this is the same fit as far as the code under test is concerned.
.statedepfit <- function(){
  utils::data('ctstantestfit', package = 'ctsem', envir = environment())
  fit <- ctstantestfit
  row <- which(fit$setup$matsetup$matrix == 3 & fit$setup$matsetup$row == 1 &
      fit$setup$matsetup$col == 1)
  fit$setup$matsetup$parname[row] <- 'log1p(exp(param)) * state[2]'
  fit
}


test_that('a linear model says nothing about state dependence, anywhere', {
  utils::data('ctstantestfit', package = 'ctsem', envir = environment())
  m <- .netmodel()

  # The point of the whole feature is that it costs a linear model nothing.
  expect_silent(net <- ctNetwork(ctstantestfit, dt = 1))
  expect_null(attr(net, 'stateDependent'))
  expect_equal(attr(net, 'stateLabel'), ctsem:::.ctContextPopLabel)
  expect_null(ctsem:::.ctNetworkStateDependentCells(ctstantestfit))
  expect_output(print(net), 'ctNetwork from a fit')
  expect_false(any(grepl('linearisation', capture.output(print(net)))))

  # And no subtitle: the plot object is the one it was before state= existed.
  g <- suppressMessages(ctNetworkPlot(ctstantestfit, dt = 1))
  expect_false('subtitle' %in% names(g$labels))
  g <- suppressMessages(ctNetworkPlot(m, dt = 1))
  expect_false('subtitle' %in% names(g$labels))
})


test_that('a state dependent fit is reported once, named, and captioned', {
  fit <- .statedepfit()
  cells <- ctsem:::.ctNetworkStateDependentCells(fit)
  expect_equal(cells$matrix, 'DRIFT')
  expect_equal(cells$kind, 'state')

  expect_message(net <- ctNetwork(fit, dt = 1), 'linearisation')
  # Named in the same words as ctSummaryMatrices, summary and ctSubjectPars,
  # which is the whole reason .ctResolveState returns a label at all.
  expect_message(ctNetwork(fit, dt = 1), ctsem:::.ctContextPopLabel, fixed = TRUE)
  # And pointed at the functions that show the variation rather than one slice.
  expect_message(ctNetwork(fit, dt = 1), 'ctPhasePortrait', fixed = TRUE)
  expect_message(ctNetwork(fit, dt = 1), 'ctStateDependencePlot', fixed = TRUE)

  expect_equal(nrow(attr(net, 'stateDependent')), 1)
  expect_equal(attr(net, 'stateLabel'), ctsem:::.ctContextPopLabel)
  expect_output(print(net), 'state dependent DRIFT')

  # Once, not once per matrix and not once per interval.
  msgs <- testthat::capture_messages(ctNetwork(fit, dt = 1))
  expect_equal(sum(grepl('linearisation', msgs)), 1)
  msgs <- testthat::capture_messages(ctNetworkPlot(fit, dt = c(.5, 2),
    networks = 'temporal'))
  expect_equal(sum(grepl('linearisation', msgs)), 1)

  expect_silent(ctNetwork(fit, dt = 1, quiet = TRUE))

  # The figure carries the caveat, because a saved plot outlives the session.
  g <- suppressMessages(ctNetworkPlot(fit, dt = 1))
  expect_true(grepl('linearisation', g$labels$subtitle))
  expect_true(grepl(ctsem:::.ctContextPopLabel, g$labels$subtitle, fixed = TRUE))
})


test_that('matsetup stateref is read, and a carrier index is not state dependence', {
  # `stateref` is the model writer's own record of a cell that materialises from
  # a state, so it is read as a second source alongside the expression scan.
  fit <- structure(list(
    standata = list(nlatent = 2L),
    setup = list(matsetup = data.frame(
      parname = c('d11', 'd22'), row = 1:2, col = 1:2, matrix = c(3L, 3L),
      stateref = c(2L, 0L)))),
    class = c('ctStanFit', 'ctFit'))
  cells <- ctsem:::.ctNetworkStaterefCells(fit)
  expect_equal(cells$matrix, 'DRIFT')
  expect_equal(cells$row, 1L)

  # A reference above nlatent is an individually varying parameter's carrier
  # state, which IS that parameter -- reporting it as state dependence would
  # tell every multilevel user their network is a linearisation when it is not.
  fit$setup$matsetup$stateref <- c(9L, 0L)
  expect_null(ctsem:::.ctNetworkStaterefCells(fit))
  expect_null(ctsem:::.ctNetworkStateDependentCells(fit))

  # A state dependent CINT is real but cannot move an edge in any of these four
  # networks, so it is not reported here.
  fit$setup$matsetup$stateref <- c(0L, 0L)
  fit$setup$matsetup$parname <- c('d11', 'c1 * state[1]')
  fit$setup$matsetup$matrix <- c(3L, 7L)
  expect_true(nrow(ctsem:::.ctFitConditionalCells(fit)) > 0)
  expect_null(ctsem:::.ctNetworkStateDependentCells(fit))
})


test_that('state= is honoured where it can be and refused where it cannot', {
  utils::data('ctstantestfit', package = 'ctsem', envir = environment())

  # The default must be indistinguishable from not passing it at all, or every
  # existing result moves.
  a <- suppressMessages(ctNetwork(ctstantestfit, dt = 1.4))
  b <- suppressMessages(ctNetwork(ctstantestfit, dt = 1.4, state = 'T0MEANS'))
  expect_identical(a$temporal, b$temporal)
  expect_identical(a$contemporaneous, b$contemporaneous)
  expect_identical(a$edges, b$edges)
  expect_identical(suppressMessages(ctSummaryMatrices(ctstantestfit)),
    suppressMessages(ctSummaryMatrices(ctstantestfit, state = 'T0MEANS')))

  # Only the julia engine can re-materialise the matrices somewhere else. This
  # used to reach ctSummaryMatrices.ctStanFit through `...` and be dropped, so
  # the same call reported a chosen state on julia and T0MEANS on stan with
  # nothing said either way.
  expect_error(ctSummaryMatrices(ctstantestfit, state = 'mean'), "backend='julia'")
  expect_error(ctSummaryMatrices(ctstantestfit, state = c(0, 0)), "backend='julia'")
  expect_error(ctNetwork(ctstantestfit, dt = 1, state = 'asymptotic'), "backend='julia'")
  expect_error(ctNetworkPlot(ctstantestfit, dt = 1, state = 'mean'), "backend='julia'")

  # state= on an already-built ctNetwork cannot be honoured by relabelling.
  net <- suppressMessages(ctNetwork(ctstantestfit, dt = 1))
  expect_message(ctNetworkPlot(net, state = 'mean'), 'state is ignored')
})


test_that('a state dependent specification says which cells it left out, in the figure too', {
  msd <- suppressMessages(suppressWarnings(ctModel(type = 'ct', n.latent = 2,
    n.manifest = 1, manifestNames = 'Y1', latentNames = c('eta1', 'eta2'),
    LAMBDA = matrix(c(1, 0), 1, 2),
    DRIFT = matrix(c('-2*log1p(exp(-2*eta2))', 0, 0, -.00001), 2, 2),
    DIFFUSION = matrix(c('diff', 0, 0, 0), 2, 2))))
  # The specification path omits such a cell rather than linearising it, so
  # there is no evaluation point to name -- the caveat is about absence.
  expect_message(ctNetwork(msd, dt = 1), 'shown as absent')
  g <- suppressMessages(ctNetworkPlot(msd, dt = 1, networks = 'temporal'))
  expect_true(grepl('absent from these edges', g$labels$subtitle))
  expect_null(attr(suppressMessages(ctNetwork(msd, dt = 1)), 'stateLabel'))
})
