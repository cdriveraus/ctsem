# ctCov(): a covariance written as the cells a model reads.
#
# The cells are a lower triangle of standard deviations and correlation
# coordinates, and what a coordinate means depends on covmattransform and on
# the other coordinates. ctCov() inverts the construction; the model remembers
# the covariance so the cells are rewritten for whichever construction reads
# them. The R forward maps it inverts are checked against the engine's own
# `sdcovsqrt2cov`, which is the construction every fit uses.

.cc_S <- matrix(c(1, .6, .3, .6, 2, .5, .3, .5, 1.5), 3,
  dimnames = list(c("a", "b", "c"), c("a", "b", "c")))

test_that("ctCov round-trips under every construction", {
  for (tf in c("z", "rawcorr", "cholesky")) {
    cells <- ctCov(.cc_S, tf)
    expect_s3_class(cells, "ctCov")
    expect_equal(c(ctsem:::.ctCovForward(unclass(cells), tf)), c(.cc_S),
      tolerance = 1e-10)
    expect_true(all(unclass(cells)[upper.tri(cells)] == 0))
    expect_equal(dimnames(cells), dimnames(.cc_S))
  }
  # Two variables under 'z': the coordinate is Fisher's z.
  two <- matrix(c(1, .5, .5, 1), 2)
  expect_equal(unclass(ctCov(two))[2, 1], atanh(.5), tolerance = 1e-10)
  expect_error(ctCov(matrix(c(1, 2, 2, 1), 2)), "not positive definite")
  # A covariance the construction does not reach is refused, not approximated.
  local_mocked_bindings(.ctCovForward = function(cells, covmattransform)
    diag(nrow(cells)), .package = "ctsem")
  expect_error(ctCov(.cc_S, "rawcorr"), "cannot be written")
})

test_that("the R constructions are the engine's", {
  skip_without_julia()
  m <- suppressMessages(ctModel(LAMBDA = diag(2), Tpoints = 3, silent = TRUE))
  spec <- suppressMessages(ctFit(data.frame(id = rep(1:2, each = 3),
    time = rep(0:2, 2), Y1 = rnorm(6), Y2 = rnorm(6)), m, backend = "julia",
    fit = FALSE))
  ctsem:::.ctJuliaModule(spec$project)
  set.seed(3)
  for (code in 0:2) {
    cells <- matrix(rnorm(16), 4)
    cells[upper.tri(cells)] <- 0
    diag(cells) <- abs(diag(cells)) + .3
    engine <- JuliaConnectoR::juliaLet(
      "Matrix(ContinuousTimeSEM.sdcovsqrt2cov(m, c))", m = cells, c = code)
    expect_equal(ctsem:::.ctCovForward(cells,
      c("rawcorr", "cholesky", "z")[code + 1L]), engine, tolerance = 1e-12)
  }
})

test_that("a ctCov matrix keeps its covariance when the construction changes", {
  m <- suppressMessages(ctModel(LAMBDA = diag(3), DIFFUSION = ctCov(.cc_S),
    Tpoints = 3, silent = TRUE))
  expect_named(m$covinput, "DIFFUSION")
  held <- function(model, tf) ctsem:::.ctCovForward(
    ctsem:::.ctCovCurrentCells(model, "DIFFUSION", 3), tf)
  expect_equal(held(m, "z"), unname(.cc_S), tolerance = 1e-10)
  # covmattransform changed after the matrix was written.
  raw <- m
  raw$covmattransform <- "rawcorr"
  raw <- ctsem:::.ctCovRefresh(raw)
  expect_equal(held(raw, "rawcorr"), unname(.cc_S), tolerance = 1e-10)
  # The R generator reads a Cholesky factor.
  expect_equal(held(ctsem:::.ctCovRefresh(m, "cholesky"), "cholesky"),
    unname(.cc_S), tolerance = 1e-10)
  # A cell edited by hand is the user's statement: kept, and the covariance
  # behind the old cells forgotten.
  edited <- m
  edited$pars$value[edited$pars$matrix %in% "DIFFUSION" &
    edited$pars$row == 2 & edited$pars$col == 1] <- 0
  edited$covmattransform <- "rawcorr"
  edited <- ctsem:::.ctCovRefresh(edited)
  expect_equal(ctsem:::.ctCovCurrentCells(edited, "DIFFUSION", 3)[2, 1], 0)
  expect_null(edited$covinput)
  # The matrix view takes one too.
  m$matrices$T0VAR <- ctCov(diag(c(1, 2, 3)))
  expect_setequal(names(m$covinput), c("DIFFUSION", "T0VAR"))
})

test_that("RAWPOPVAR <- ctCov() is the Laplace population covariance", {
  skip_without_julia()
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-.5), DIFFUSION = matrix(.5), MANIFESTVAR = matrix(.3),
    MANIFESTMEANS = matrix("mm"), T0VAR = matrix(1), T0MEANS = matrix(0),
    CINT = matrix("cint"))))
  m$pars$indvarying <- m$pars$param %in% c("mm", "cint")
  P <- matrix(c(.5, .2, .2, .3), 2, dimnames = list(c("mm", "cint"), c("mm", "cint")))
  m$RAWPOPVAR <- ctCov(P)
  expect_error(m$RAWPOPVAR <- ctCov(matrix(1, 1, 1, dimnames = list("mm", "mm"))),
    "whole matrix")
  spec <- suppressMessages(ctFit(data.frame(id = rep(1:6, each = 3),
    time = rep(0:2, 6), Y1 = rnorm(18)), m, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  popcov <- as.matrix(ctsem:::.ctBackendJuliaValue(
    ctsem:::.ctJuliaModule(spec$project)$ctsem_laplace_popcov(
      ctsem:::.ctJuliaObjective(handle),
      ctsem:::.ctJuliaNumericVector(numeric(max(1L, ctsem:::.ctBackendNpar(spec)))),
      1L)))
  order <- match(spec$laplace$levels[[1]]$param, rownames(P))
  expect_equal(popcov, unname(P[order, order]), tolerance = 1e-12,
    ignore_attr = TRUE)
})
