# Shared dimensions in blocks at a reduced-rank level: ctFit(popblocks =)
# (R/ctPopBlocks.R; `load_row`, `load_col` in the engine).
#
# Two independent answers come first and fit nothing. One block naming every
# effect is the triangle a level without blocks has, so its objective is the
# plain poprank model's to rounding: that checks the pattern's plumbing. And
# effects in no block, with own variances, are the same model as those
# effects on a second id level identical to the subject at rank 0: a
# different route through the engine to the same integral. Then the pattern's
# point, exact zero covariance between effects with no block in common, and
# one fit with what reports it.

.blocks_data <- function(nsub = 40L, nt = 10L, seed = 5L) {
  set.seed(seed)
  lam <- c(1, 0.8, 1.2)
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    times <- cumsum(c(0, stats::runif(nt - 1L, 0.5, 1.5)))
    zl <- stats::rnorm(1); zd <- stats::rnorm(1); zi <- stats::rnorm(1)
    a <- -0.6 + 0.15 * zd
    cint <- 0.6 * zl - 0.1 * zd
    mm <- c(0, 0.4 * zi + 0.2 * zd, -0.3 * zi + stats::rnorm(1, 0, 0.2))
    eta <- numeric(nt); eta[1] <- stats::rnorm(1)
    for (t in 2:nt) {
      e <- exp(a * (times[t] - times[t - 1]))
      eta[t] <- e * eta[t - 1] + (e - 1) / a * cint +
        stats::rnorm(1, 0, sqrt((1 - e^2) / (-2 * a)))
    }
    y <- sapply(1:3, function(j) lam[j] * eta + mm[j] + stats::rnorm(nt, 0, 0.3))
    data.frame(id = i, id2 = i, time = times, Y1 = y[, 1], Y2 = y[, 2], Y3 = y[, 3])
  }))
}

.blocks_model <- function(ids = "id", lev = "lev||TRUE", dr = "dr||TRUE",
  mm2 = "mm2||TRUE", mm3 = "mm3||TRUE") {
  suppressMessages(ctModel(type = "ct", manifestNames = c("Y1", "Y2", "Y3"),
    latentNames = "eta", id = ids, LAMBDA = matrix(c(1, 0.8, 1.2)),
    DRIFT = matrix(dr), CINT = matrix(lev), DIFFUSION = matrix("q"),
    MANIFESTVAR = matrix(c("e1", 0, 0, 0, "e2", 0, 0, 0, "e3"), 3),
    MANIFESTMEANS = matrix(c(0, mm2, mm3)), T0MEANS = matrix(0), T0VAR = matrix(1)))
}

.blocks2 <- list(level = list(rank = 1, effects = c("lev", "dr")),
  items = list(rank = 1, effects = c("mm2", "mm3", "dr")))

.blocks_prep <- function(m, data, ...) suppressMessages(ctFit(data, m,
  backend = "julia", intoverpop = "laplace", fit = FALSE, cores = 2, ...))

test_that("blocks are read, laid out and refused where they cannot apply", {
  m <- .blocks_model()
  levelcolumns <- c(id = "indvarying")
  b <- ctsem:::.ctPopblocksResolve(m, .blocks2, "id", levelcolumns)
  expect_equal(names(b), "id")
  expect_equal(vapply(b$id, `[[`, character(1), "name"), c("level", "items"))
  expect_setequal(b$id[[2]]$param, c("mm2", "mm3", "dr"))
  # Block `b` takes the next `rank` dimensions, its effects in level order
  # with the first `rank` its basis: one block over every effect is the
  # triangle.
  names <- c("dr", "mm2", "mm3", "lev")
  pattern <- ctsem:::.ctJuliaBlockPattern(b$id, names)
  expect_equal(pattern$row, c(1L, 4L, 1L, 2L, 3L))
  expect_equal(pattern$col, c(1L, 1L, 2L, 2L, 2L))
  expect_equal(pattern$prior, 1 / sqrt(c(2, 1, 2, 1, 1)))
  one <- ctsem:::.ctJuliaBlockPattern(list(list(name = "all", rank = 2L,
    param = names)), names)
  tri <- ctsem:::.ctJuliaLoadingEntries(list(nrandom = 4L, rank = 2L,
    load_index = 1:7))
  expect_equal(one$row, tri$row)
  expect_equal(one$col, tri$col)

  expect_error(ctsem:::.ctPopblocksResolve(m, list(a = list(rank = 2,
    effects = "lev")), "id", levelcolumns), "at most their number")
  expect_error(ctsem:::.ctPopblocksResolve(m, list(a = list(rank = 1,
    effects = c("lev", "dr")), b = list(rank = 1, effects = c("dr", "lev"))),
    "id", levelcolumns), "same effects")
  expect_error(ctsem:::.ctPopblocksResolve(m, list(a = list(rank = 1,
    effects = "nope")), "id", levelcolumns), "nope")
  data <- .blocks_data(nsub = 6L, nt = 3L)
  expect_error(.blocks_prep(m, data, popblocks = list(level = list(rank = 1,
    effects = "lev"))), "in no block and without an own variance")
  expect_error(.blocks_prep(m, data, popblocks = .blocks2, poprank = 3),
    "its blocks give its rank")
  # `'auto'`, the default written out, and NULL, the argument left out, defer
  # to the blocks.
  for (rank in list("auto", NULL)) {
    p <- .blocks_prep(m, data, popblocks = .blocks2, poprank = rank)
    expect_equal(p$laplace$levels[[1]]$rank, 2L)
  }
  expect_error(suppressMessages(ctFit(data, m, backend = "julia", fit = FALSE,
    intoverpop = "augmented", popblocks = .blocks2)), "intoverpop='laplace'")
})

test_that("one block is the triangle, and no block is a level of its own", {
  skip_without_julia()
  data <- .blocks_data()
  p0 <- .blocks_prep(.blocks_model(), data, poprank = 2)
  p1 <- .blocks_prep(.blocks_model(), data,
    popblocks = list(all = list(rank = 2, effects = TRUE)))
  expect_equal(ctsem:::.ctBackendNpar(p1), ctsem:::.ctBackendNpar(p0))
  set.seed(1)
  theta <- stats::rnorm(ctsem:::.ctBackendNpar(p0), 0, 0.4)
  e0 <- ctJuliaEvaluate(p0, theta, cores = 2)
  e1 <- ctJuliaEvaluate(p1, theta, cores = 2)
  expect_equal(e1$value, e0$value, tolerance = 1e-12)
  expect_equal(e1$gradient, e0$gradient, tolerance = 1e-10)

  # Effects in no block with own variances, against the same effects on a
  # second id level identical to the subject, at rank 0.
  pa <- .blocks_prep(.blocks_model(), data, popblocks = list(level = list(
    rank = 1, effects = c("lev", "dr"))), popunique = c("mm2", "mm3"))
  pb <- .blocks_prep(.blocks_model(ids = c("id", "id2"), lev = "lev||id",
    dr = "dr||id", mm2 = "mm2||id2", mm3 = "mm3||id2"), data,
    poprank = c(id = 1, id2 = 0), popunique = c(id2 = TRUE))
  na <- ctsem:::.ctBackendRawParameterNames(list(model_spec = pa),
    ctsem:::.ctBackendNpar(pa))
  nb <- sub("\\.(id|id2)$", "", ctsem:::.ctBackendRawParameterNames(
    list(model_spec = pb), ctsem:::.ctBackendNpar(pb)))
  expect_setequal(na, nb)
  set.seed(3)
  theta <- stats::rnorm(length(na), 0, 0.4)
  ea <- ctJuliaEvaluate(pa, theta, cores = 2)
  eb <- ctJuliaEvaluate(pb, theta[match(nb, na)], cores = 2)
  expect_equal(ea$value, eb$value, tolerance = 1e-8)
  expect_equal(ea$gradient[match(nb, na)], eb$gradient, tolerance = 1e-6)
})

test_that("effects with no block in common are uncorrelated, fitted and reported", {
  skip_without_julia()
  data <- .blocks_data()
  p <- .blocks_prep(.blocks_model(), data, popblocks = .blocks2)
  lv <- p$laplace$levels[[1]]
  set.seed(2)
  theta <- stats::rnorm(ctsem:::.ctBackendNpar(p), 0, 0.5)
  module <- ctsem:::.ctJuliaModule(p$project)
  S <- as.matrix(module$ctsem_laplace_popcov(ctsem:::.ctJuliaObjective(p,
    hybrid = TRUE), ctsem:::.ctJuliaNumericVector(theta), 1L))
  dimnames(S) <- list(lv$param, lv$param)
  expect_equal(unname(S["lev", c("mm2", "mm3")]), c(0, 0))
  expect_true(all(S["dr", c("lev", "mm2", "mm3")] != 0))

  fit <- suppressMessages(ctFit(data, .blocks_model(), backend = "julia",
    intoverpop = "laplace", cores = 2, popblocks = .blocks2))
  s <- summary(fit)
  expect_equal(unname(s$rawpopcorr[c("lev__mm2", "lev__mm3"), "mean"]), c(0, 0))
  expect_match(s$popsdNote, "come in blocks", fixed = TRUE)
  rawnames <- ctsem:::.ctBackendRawParameterNames(fit, length(fit$estimate$raw))
  expect_true(all(c("poploading_lev_dim1", "poploading_mm2_dim2") %in% rawnames))
  expect_false("poploading_lev_dim2" %in% rawnames)
  expect_equal(dim(ctSubjectPars(fit))[2L], length(unique(data$id)))
  # SAEM's loading moves keep the pattern: its phase ends where the fit does.
  saem <- suppressMessages(ctFit(data, .blocks_model(), backend = "julia",
    intoverpop = "laplace", cores = 2, popblocks = .blocks2,
    optimcontrol = list(saem = TRUE)))
  expect_equal(saem$estimate$loglik, fit$estimate$loglik, tolerance = 1e-3)
})
