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
  local_mocked_bindings(.ctCovForward = function(cells, covmattransform, zero = NULL)
    diag(nrow(cells)), .package = "ctsem")
  expect_error(ctCov(.cc_S, "rawcorr"), "cannot be written")
})

# Under 'z' a cell of T0VAR, DIFFUSION or MANIFESTVAR fixed at zero is a zero
# covariance. Between variables that other cells connect the construction
# solves the cell to make it one.
.cc_pattern <- function(d, share, sdcoord) {
  zero <- matrix(FALSE, d, d)
  zero[lower.tri(zero)] <- stats::runif(d * (d - 1) / 2) < share
  cells <- diag(exp(stats::rnorm(d, 0, .3)), d)
  cells[lower.tri(cells)] <- stats::rnorm(d * (d - 1) / 2, 0, sdcoord)
  cells[zero] <- 0
  list(cells = cells, zero = zero)
}

test_that("under 'z' a cell of zero is a zero covariance", {
  forward <- ctsem:::.ctCovForward
  cells <- diag(3)
  cells[3, 1:2] <- 0.5
  zero <- ctsem:::.ctCovZeroCells(cells)
  expect_identical(which(zero), 2L)
  S <- forward(cells, "z", zero)
  expect_identical(S[2, 1], 0)
  expect_equal(S[3, 1:2], rep(0.4306862, 2), tolerance = 1e-6)
  # Held at 0 as a coordinate it is not one, and the other constructions read
  # the cell as it is.
  expect_equal(forward(cells, "z")[2, 1], 0.1152759, tolerance = 1e-6)
  expect_identical(forward(cells, "rawcorr", zero), forward(cells, "rawcorr"))
  # Between separate blocks nothing is solved.
  blocks <- diag(4)
  blocks[2, 1] <- 0.7
  blocks[4, 3] <- -0.4
  expect_false(any(ctsem:::.ctCovLinkedZeros(ctsem:::.ctCovZeroCells(blocks))))
  expect_equal(forward(blocks, "z", ctsem:::.ctCovZeroCells(blocks)),
    forward(blocks, "z"), tolerance = 1e-14)
  # Any pattern: exact zeros, positive definite, and the same matrix whatever
  # the order of the variables.
  set.seed(4)
  for (rep in 1:25) {
    p <- .cc_pattern(6, 0.4, 1)
    S <- forward(p$cells, "z", p$zero)
    expect_true(all(S[p$zero] == 0))
    expect_gt(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values), 0)
    perm <- sample(6)
    full <- p$cells + t(p$cells) - diag(diag(p$cells))
    cellsp <- full[perm, perm]
    cellsp[upper.tri(cellsp)] <- 0
    zerop <- (p$zero | t(p$zero))[perm, perm] & lower.tri(p$zero)
    expect_equal(forward(cellsp, "z", zerop), S[perm, perm], tolerance = 1e-10)
  }
  # Thirty variables.
  p <- .cc_pattern(30, 0.3, 0.3)
  S <- forward(p$cells, "z", p$zero)
  expect_true(all(S[p$zero] == 0))
  expect_gt(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("ctCov writes a zero covariance as a cell of zero", {
  forward <- ctsem:::.ctCovForward
  S <- matrix(c(1, 0, .4, 0, 1, .4, .4, .4, 1), 3)
  cells <- unclass(ctCov(S))
  attributes(cells) <- list(dim = dim(cells))
  expect_identical(cells[2, 1], 0)
  expect_equal(forward(cells, "z", ctsem:::.ctCovZeroCells(cells)), S, tolerance = 1e-10)
  # A covariance that is not zero is never written as a cell of exactly zero,
  # which would be read as one.
  coord <- diag(3)
  coord[3, 1:2] <- 0.5
  S2 <- forward(coord, "z")
  back <- ctsem:::.ctCovCells(S2, "z")
  expect_false(back[2, 1] == 0)
  expect_equal(forward(back, "z", ctsem:::.ctCovZeroCells(back)), S2, tolerance = 1e-10)
  # The population's cells are all coordinates.
  pop <- ctsem:::.ctCovCells(S, "z", zeros = FALSE)
  expect_lt(pop[2, 1], -0.05)
  expect_equal(forward(pop, "z"), S, tolerance = 1e-10)
  # Random targets with a pattern of zeros are reached.
  set.seed(9)
  for (rep in 1:10) {
    p <- .cc_pattern(5, 0.4, 0.6)
    target <- forward(p$cells, "z", p$zero)
    got <- ctsem:::.ctCovCells(target, "z")
    expect_equal(forward(got, "z", ctsem:::.ctCovZeroCells(got)), target,
      tolerance = 1e-9)
  }
})

test_that("a model's cells fixed at zero are found, and stan refuses them", {
  m <- suppressMessages(ctModel(LAMBDA = diag(3), Tpoints = 5, silent = TRUE,
    DIFFUSION = matrix(c("d11", 0, "d31", 0, "d22", "d32", 0, 0, "d33"), 3)))
  expect_equal(unname(which(ctsem:::.ctCovSolvedCells(m, "DIFFUSION"), arr.ind = TRUE)),
    matrix(2:1, 1))
  # the default MANIFESTVAR is diagonal: zeros between separate variables
  expect_false(any(ctsem:::.ctCovSolvedCells(m, "MANIFESTVAR")))
  expect_error(ctsem:::.ctCovRefuseSolved(m), "DIFFUSION\\[2,1\\].*backend = 'julia'")
  raw <- m
  raw$covmattransform <- "rawcorr"
  expect_false(any(ctsem:::.ctCovSolvedCells(raw, "DIFFUSION")))
  expect_silent(ctsem:::.ctCovRefuseSolved(raw))
  data <- data.frame(id = rep(1:3, each = 4), time = rep(0:3, 3), Y1 = rnorm(12),
    Y2 = rnorm(12), Y3 = rnorm(12))
  expect_error(suppressMessages(ctFit(data, m, backend = "stan", fit = FALSE)),
    "DIFFUSION\\[2,1\\]")
})

.cc_zero_model <- function() suppressMessages(ctModel(LAMBDA = diag(3), silent = TRUE,
  manifestNames = paste0("Y", 1:3), latentNames = paste0("eta", 1:3),
  DRIFT = matrix(c("dr1", 0, 0, 0, "dr2", 0, "cr31", 0, "dr3"), 3, byrow = TRUE),
  DIFFUSION = matrix(c("d11", 0, "d31", 0, "d22", "d32", 0, 0, "d33"), 3),
  T0VAR = matrix(c("t11", "t21", 0, 0, "t22", "t32", 0, 0, "t33"), 3),
  MANIFESTVAR = matrix(c("m11", "m21", "m31", 0, "m22", 0, 0, 0, "m33"), 3),
  T0MEANS = matrix(0, 3, 1), MANIFESTMEANS = matrix(0, 3, 1), CINT = matrix(0, 3, 1)))

.cc_zero_data <- function(nsub = 20L, nobs = 8L) {
  set.seed(21)
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    x <- matrix(0, nobs, 3)
    x[1, ] <- stats::rnorm(3)
    for (t in 2:nobs) {
      shared <- stats::rnorm(1)
      x[t, ] <- 0.6 * x[t - 1, ] + c(0.5 * shared + stats::rnorm(1, 0, .6),
        -0.4 * shared + stats::rnorm(1, 0, .6), shared)
    }
    data.frame(id = i, time = seq_len(nobs) - 1, Y1 = x[, 1] + stats::rnorm(nobs, 0, .3),
      Y2 = x[, 2] + stats::rnorm(nobs, 0, .3), Y3 = x[, 3] + stats::rnorm(nobs, 0, .3))
  }))
}

test_that("the engine solves the cells R does", {
  skip_without_julia()
  spec <- suppressMessages(ctFit(.cc_zero_data(4L, 4L), .cc_zero_model(),
    backend = "julia", fit = FALSE))
  ctsem:::.ctJuliaModule(spec$project)
  construct <- JuliaConnectoR::juliaEval(paste0("(m, r, c) -> Matrix(",
    "ContinuousTimeSEM.sdcovsqrt2cov(m, ContinuousTimeSEM.CTSEMCovCode(2, ",
    "ContinuousTimeSEM.CTSEMZeroCells(Int.(vec(r)), Int.(vec(c))))))"))
  set.seed(6)
  for (rep in 1:8) {
    p <- .cc_pattern(if (rep == 8) 30 else 6, 0.4, if (rep == 8) 0.3 else 0.8)
    linked <- which(ctsem:::.ctCovLinkedZeros(p$zero), arr.ind = TRUE)
    if (nrow(linked) < 2L) next
    engine <- construct(p$cells, matrix(linked[, 1]), matrix(linked[, 2]))
    expect_equal(ctsem:::.ctCovForward(p$cells, "z", p$zero), engine, tolerance = 1e-10)
  }
})

test_that("a model with zero covariances is fitted with them", {
  skip_without_julia()
  model <- .cc_zero_model()
  data <- .cc_zero_data()
  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
  set.seed(8)
  raw <- stats::rnorm(max(spec$parameter_table$parnumber, na.rm = TRUE), 0, .4)
  at <- ctJuliaEvaluate(spec, raw, gradient = TRUE)
  expect_true(is.finite(at$value))
  # The reverse pass through the solved cells, against central differences of
  # the likelihood.
  finite <- vapply(seq_along(raw), function(i) {
    step <- rep(0, length(raw))
    step[i] <- 1e-5
    (ctJuliaEvaluate(spec, raw + step, gradient = FALSE)$value -
        ctJuliaEvaluate(spec, raw - step, gradient = FALSE)$value) / 2e-5
  }, numeric(1))
  expect_equal(as.numeric(at$gradient), finite, tolerance = 1e-6)
  # The matrices the engine builds from the model's table are the ones R's
  # forward map gives from the model's cells: the two find the same cells to
  # solve, and solve them alike.
  pop <- ctsem:::.ctBackendPopArrays(structure(list(model_spec = spec, model = model,
    backend = "julia", estimate = list(raw = raw, loglik = NA_real_)),
    class = c("ctJuliaFit", "ctFit")))
  for (pair in list(c("DIFFUSION", "DIFFUSIONcov"), c("MANIFESTVAR", "MANIFESTcov"),
      c("T0VAR", "T0cov"))) {
    cells <- drop(pop[[paste0("pop_", pair[1])]])
    cov <- drop(pop[[paste0("pop_", pair[2])]])
    solved <- ctsem:::.ctCovSolvedCells(model, pair[1])
    expect_identical(sum(solved), 1L, info = pair[1])
    expect_true(all(cov[solved | t(solved)] == 0), info = pair[1])
    expect_equal(unname(cov), ctsem:::.ctCovForward(unname(cells), "z", solved),
      tolerance = 1e-10, info = pair[1])
  }
  # What a fit reports: exact zeros where the model has them, covariances
  # elsewhere.
  fit <- suppressWarnings(suppressMessages(ctFit(data, model, backend = "julia",
    cores = 2)))
  mats <- ctSummaryMatrices(fit)
  expect_identical(unname(mats$DIFFUSIONcov[2, 1]), 0)
  expect_identical(unname(mats$T0cov[3, 1]), 0)
  expect_identical(unname(mats$MANIFESTcov[3, 2]), 0)
  expect_gt(abs(mats$DIFFUSIONcov[3, 1]), 0.01)
  expect_gt(abs(mats$DIFFUSIONcov[3, 2]), 0.01)
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
