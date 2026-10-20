# ctHalfLife() against closed forms: an Ornstein-Uhlenbeck process halves its
# autocorrelation at log(2) / rate, and a discrete one with coefficient a at
# log(0.5) / log(a) occasions.

test_that("half-lives match the closed forms", {
  one <- ctsem:::.ctHalfLifeOne(matrix(-0.5), matrix(1), 1L, NULL, TRUE)
  expect_equal(one, log(2) / 0.5, tolerance = 1e-8)
  discrete <- ctsem:::.ctHalfLifeOne(matrix(0.8), matrix(1), 1L, NULL, FALSE)
  expect_equal(discrete, log(0.5) / log(0.8), tolerance = 1e-8)
  # A fast process pushed around by a slow one: its whole autocorrelation
  # holds on far longer than its own part does, whose half-life is its own
  # rate's.
  A <- matrix(c(-1, 0, 0.5, -0.01), 2, 2)
  G <- diag(c(1, 0.1))
  whole <- ctsem:::.ctHalfLifeOne(A, G, 1L, NULL, TRUE)
  own <- ctsem:::.ctHalfLifeOne(A, G, 1L, 1L, TRUE)
  expect_equal(own, log(2), tolerance = 1e-8)
  expect_gt(whole, 5 * own)
  # The slow process is not reached by the fast one's innovations.
  expect_true(is.na(ctsem:::.ctHalfLifeOne(A, G, 2L, 1L, TRUE)))
  # A source whose innovations correlate with the others' has no part of its own.
  expect_error(ctsem:::.ctHalfLifeOne(A, matrix(c(1, .05, .05, .1), 2), 1L, 1L,
    TRUE), "correlated")
  # No stationary distribution, no half-life.
  expect_true(is.na(ctsem:::.ctHalfLifeOne(matrix(0.1), matrix(1), 1L, NULL, TRUE)))
})

test_that("ctHalfLife gives a draw per sample and process", {
  skip_on_cran()
  hl <- suppressMessages(ctHalfLife(ctstantestfit, nsamples = 20))
  expect_equal(dim(hl)[c(2L, 3L)], c(1L, 2L))
  expect_equal(dimnames(hl)$Process, ctstantestfit$ctstanmodel$latentNames)
  expect_true(all(is.finite(hl) & hl > 0))
  expect_error(ctHalfLife(ctstantestfit, processes = 'nope'), "does not have")
})
