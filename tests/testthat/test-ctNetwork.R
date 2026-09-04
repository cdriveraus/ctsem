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
