# The bound on random-effect correlations (R/ctCorBound.R, the engine's
# `_ctsem_corbound`): each effect's multiple correlation with the others at its
# level at most 0.99, held by a wall in the objective on every julia route.
#
# The independent answer is the wall computed here, in R, from the population
# correlation matrix `.ctCovForward` builds: the objective with the bound minus
# the objective without it. Past the bound the two routes must then agree,
# which they did not while only the laplace route clamped the coordinates, and
# the adjoint gradient must agree with forward mode. One fit, on data whose
# maximum-likelihood population covariance is singular (rank 2 of 4), must end
# on the bound, say so, and give every correlation a width.

.cb_model <- function() ctModel(type = "ct", latentNames = c("eta1", "eta2"),
  manifestNames = c("Y1", "Y2"), MANIFESTVAR = diag(.1, 2),
  TDpredNames = "TD1", TIpredNames = c("TI1", "TI2", "TI3"), LAMBDA = diag(2),
  silent = TRUE)

.cb_spec <- function(route) structure(suppressWarnings(suppressMessages(ctFit(
  ctstantestdat, .cb_model(), backend = "julia", intoverpop = route, priors = FALSE,
  fit = FALSE))), class = c("ctJuliaModel", "ctFitModel"))

# The wall, from the correlation matrix the block's coordinates build.
.cb_wall <- function(block, raw, bound = 0.99, steepness = 100) {
  cells <- diag(block$k)
  low <- which(lower.tri(cells), arr.ind = TRUE)
  low <- low[order(low[, 2L], low[, 1L]), , drop = FALSE]
  cells[low] <- ifelse(block$index > 0L, raw[pmax(block$index, 1L)], block$fixed)
  R <- ctsem:::.ctCovForward(cells, if (block$code == 2L) "z" else "rawcorr")
  excess <- log(diag(R) * diag(solve(R))) + log1p(-bound^2)
  -steepness * sum(pmax(excess, 0)^3)
}

test_that("the bound is the same wall on both routes, and its gradient is right", {
  skip_without_julia()
  specs <- list(augmented = .cb_spec("augmented"), laplace = .cb_spec("laplace"))
  blocks <- lapply(specs, ctsem:::.ctCorBoundSpec)
  expect_length(blocks$augmented, 1L)
  expect_identical(blocks$augmented[[1]][c("k", "code", "index")],
    blocks$laplace[[1]][c("k", "code", "index")])
  block <- blocks$laplace[[1]]
  names <- lapply(specs, function(s) as.character(ctsem:::.ctBackendRawParameterNames(
    list(model_spec = s), ctsem:::.ctBackendNpar(s))))
  set.seed(4)
  base <- stats::setNames(stats::rnorm(length(names$laplace), 0, .2), names$laplace)
  # Inside the bound, and past it as far as an optimiser goes: 1.2 times the
  # coordinates a bounded fit of these data ends at, where every multiple
  # correlation is about 0.995. Much further out the correlation matrix is
  # near enough to singular that the two routes' likelihoods agree only to
  # about 1e-7.
  inside <- base
  past <- base
  past[block$index] <- 1.2 * c(0.65, 0.927, 2.13, -2.02, 1.5, -0.178)
  module <- ctsem:::.ctJuliaModule(specs$laplace$project)
  on.exit(ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_set_correlation_bound!(0.99)"),
    add = TRUE)
  value <- function(s, x, route) ctJuliaEvaluate(s, unname(x[names[[route]]]),
    gradient = FALSE, cores = 1)$value
  for (x in list(inside, past)) {
    with_bound <- vapply(names(specs), function(r) value(specs[[r]], x, r), numeric(1))
    # The routes' likelihoods are the same function here (linear effects), to
    # the Laplace inner solve's tolerance at a point this far from any optimum.
    expect_equal(with_bound[["augmented"]], with_bound[["laplace"]], tolerance = 1e-6)
    ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_set_correlation_bound!(0)")
    without <- value(specs$laplace, x, "laplace")
    ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_set_correlation_bound!(0.99)")
    expect_equal(with_bound[["laplace"]] - without, .cb_wall(block, unname(x[names$laplace])),
      tolerance = 1e-8)
  }
  expect_lt(.cb_wall(block, unname(past[names$laplace])), -1)
  expect_equal(.cb_wall(block, unname(inside[names$laplace])), 0)
  for (r in names(specs)) {
    x <- unname(past[names[[r]]])
    a <- ctJuliaEvaluate(specs[[r]], x, gradient = TRUE, cores = 1)$gradient
    f <- ctJuliaEvaluate(specs[[r]], x, gradient = TRUE, gradient_method = "forward",
      cores = 1)$gradient
    expect_equal(as.numeric(a), as.numeric(f), tolerance = 1e-8)
  }
})

test_that("a fit whose covariance the data put on the singular boundary ends on the bound", {
  skip_without_julia()
  warnings <- character()
  fit <- withCallingHandlers(suppressMessages(ctFit(ctstantestdat, .cb_model(),
    backend = "julia", intoverpop = "augmented", priors = FALSE, cores = 1, verbose = 0)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("at their bound", warnings)))
  expect_false(any(grepl("do not identify", warnings)))
  expect_true(any(fit$corbound$at_bound))
  expect_equal(fit$corbound$rank[1], 2)
  raw <- fit$estimate$raw
  block <- ctsem:::.ctCorBoundSpec(fit$model_spec)[[1]]
  expect_lt(max(abs(raw[block$index])), 4)
  s <- summary(fit)
  expect_true(all(is.finite(s$rawpopcorr$sd)))
})
