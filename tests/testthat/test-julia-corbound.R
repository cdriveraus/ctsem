# The bound on random-effect correlations (R/ctCorBound.R, the engine's
# `_ctsem_corbound`): each effect's multiple correlation with the others at its
# level at most the engine's bound (`_CTSEM_CORBOUND`, 0.999), held by a wall in
# the objective on every julia route.
#
# The independent answer is the wall computed here, in R, from the population
# correlation matrix `.ctCovForward` builds: the objective minus the subjects'
# own log likelihoods, with no prior to account for. Past the bound the two routes must then agree,
# which they did not while only the laplace route clamped the coordinates, and
# the adjoint gradient must agree with forward mode. One fit, on data whose
# maximum-likelihood population covariance is singular (rank 2 of 4), must end
# on the bound, say so, and give every correlation a width.
#
# The model's own covariance matrices (T0VAR, DIFFUSION, MANIFESTVAR) are held
# by the same wall since 2026-11-08: one more block per matrix with a free
# correlation, after the population's.

.cb_model <- function() ctModel(type = "ct", latentNames = c("eta1", "eta2"),
  manifestNames = c("Y1", "Y2"), MANIFESTVAR = diag(.1, 2),
  TDpredNames = "TD1", TIpredNames = c("TI1", "TI2", "TI3"), LAMBDA = diag(2),
  silent = TRUE)

.cb_spec <- function(route) structure(suppressWarnings(suppressMessages(ctFit(
  ctstantestdat, .cb_model(), backend = "julia", intoverpop = route, priors = FALSE,
  fit = FALSE))), class = c("ctJuliaModel", "ctFitModel"))

# The wall, from the correlation matrix the block's coordinates build.
.cb_wall <- function(block, raw, bound = ctsem:::.ctCorBoundValue(), steepness = 100) {
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
  # The population's block, then DIFFUSION's: T0VAR's correlation is fixed
  # where T0MEANS vary, and MANIFESTVAR is diagonal here.
  expect_length(blocks$augmented, 2L)
  expect_null(blocks$augmented[[1]]$matrix)
  expect_identical(blocks$augmented[[2]]$matrix, "DIFFUSION")
  for (i in 1:2) expect_identical(blocks$augmented[[i]][c("k", "code", "index")],
    blocks$laplace[[i]][c("k", "code", "index")])
  block <- blocks$laplace[[1]]
  diffusion <- blocks$laplace[[2]]
  expect_equal(diffusion$k, 2L)
  names <- lapply(specs, function(s) as.character(ctsem:::.ctBackendRawParameterNames(
    list(model_spec = s), ctsem:::.ctBackendNpar(s))))
  set.seed(4)
  base <- stats::setNames(stats::rnorm(length(names$laplace), 0, .2), names$laplace)
  # Inside the bound, and past it: along the coordinates a bounded fit of
  # these data ends at, the first multiple of them whose wall is below -1, a
  # little past the bound as an optimiser goes. Much further out the
  # correlation matrix is near enough to singular that the two routes'
  # likelihoods agree only to about 1e-7.
  inside <- base
  past <- base
  direction <- c(0.65, 0.927, 2.13, -2.02, 1.5, -0.178)
  for (scale in seq(1.2, 4, by = 0.1)) {
    past[block$index] <- scale * direction
    if (.cb_wall(block, unname(past[names$laplace])) < -1) break
  }
  value <- function(s, x, route) ctJuliaEvaluate(s, unname(x[names[[route]]]),
    gradient = FALSE, cores = 1)$value
  # And DIFFUSION's correlation past its bound too: its own wall adds.
  pastdiff <- past
  for (v in seq(2, 8, by = 0.1)) {
    pastdiff[diffusion$index] <- v
    if (.cb_wall(diffusion, unname(pastdiff[names$laplace])) < -1) break
  }
  expect_lt(.cb_wall(diffusion, unname(pastdiff[names$laplace])), -1)
  for (x in list(inside, past, pastdiff)) {
    with_bound <- vapply(names(specs), function(r) value(specs[[r]], x, r), numeric(1))
    # The routes' likelihoods are the same function here (linear effects), to
    # the Laplace inner solve's tolerance at a point this far from any optimum.
    # Not where the diffusion is all but singular: there they differ by 0.02,
    # so each route's wall is checked against the one computed here instead.
    if (!identical(x, pastdiff)) {
      expect_equal(with_bound[["augmented"]], with_bound[["laplace"]], tolerance = 1e-6)
    }
    wall <- .cb_wall(block, unname(x[names$laplace])) +
      .cb_wall(diffusion, unname(x[names$laplace]))
    for (r in names(specs)) {
      without <- sum(ctJuliaEvaluate(specs[[r]], unname(x[names[[r]]]),
        gradient = FALSE, contributions = TRUE, cores = 1)$subject_loglik)
      expect_equal(with_bound[[r]] - without, wall, tolerance = 1e-8)
    }
  }
  expect_lt(.cb_wall(block, unname(past[names$laplace])), -1)
  expect_equal(.cb_wall(block, unname(inside[names$laplace])), 0)
  for (r in names(specs)) for (x in list(unname(past[names[[r]]]),
      unname(pastdiff[names[[r]]]))) {
    a <- ctJuliaEvaluate(specs[[r]], x, gradient = TRUE, cores = 1)$gradient
    f <- ctJuliaEvaluate(specs[[r]], x, gradient = TRUE, gradient_method = "forward",
      cores = 1)$gradient
    expect_equal(as.numeric(a), as.numeric(f), tolerance = 1e-8)
  }
})

test_that("the wall is finite however far out a 'z' correlation is", {
  skip_without_julia()
  # Past a coordinate of about 19, 'z' rounds a correlation to one. The wall
  # took that matrix's inverse, so the objective came back NaN at a point
  # the filter evaluates. Two effects' wall has a closed form there:
  # log((R^-1)_jj R_jj) = 2 log cosh z, for each of the two.
  spec <- .cb_spec("augmented")
  diffusion <- ctsem:::.ctCorBoundSpec(spec)[[2]]
  expect_identical(diffusion$matrix, "DIFFUSION")
  expect_equal(diffusion$code, 2L)
  set.seed(4)
  x <- stats::rnorm(ctsem:::.ctBackendNpar(spec), 0, .2)
  logcosh <- function(z) abs(z) + log1p(exp(-2 * abs(z))) - log(2)
  limit <- -log1p(-ctsem:::.ctCorBoundValue()^2)
  for (z in c(-40, 320)) {
    x[diffusion$index] <- z
    v <- ctJuliaEvaluate(spec, x, gradient = TRUE, contributions = TRUE, cores = 1)
    excess <- 2 * logcosh(z) - limit
    expect_true(all(is.finite(v$gradient)))
    expect_equal(v$value - sum(v$subject_loglik), -2 * 100 * excess^3, tolerance = 1e-10)
    # The likelihood's own slope there is sech(z)^2 times something modest.
    expect_equal(v$gradient[diffusion$index], -2 * 300 * excess^2 * 2 * tanh(z),
      tolerance = 1e-8)
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
  # A sample's draws at the estimate point to the rank the fit's warning
  # states there, which is what sampling's message suggests as poprank.
  drawn <- fit
  drawn$estimate$rawposterior <- matrix(as.numeric(raw), 3L, length(raw), byrow = TRUE)
  ranks <- ctsem:::.ctCorBoundDrawRank(drawn)
  expect_equal(ranks$rank[1], fit$corbound$rank[1])
  expect_equal(ranks$k[1], fit$corbound$k[1])
  expect_null(ctsem:::.ctCorBoundDrawRank(fit[setdiff(names(fit), "estimate")]))
  block <- ctsem:::.ctCorBoundSpec(fit$model_spec)[[1]]
  expect_lt(max(abs(raw[block$index])), 4)
  s <- summary(fit)
  expect_true(all(is.finite(s$rawpopcorr$sd)))
})

test_that("a fit whose processes share one noise source ends on DIFFUSION's bound", {
  skip_without_julia()
  # Two processes with the same drift driven by the same shocks: the
  # maximum-likelihood diffusion has correlation one.
  set.seed(11)
  a <- exp(-0.25)
  shock <- sqrt(1 - a^2)
  dat <- do.call(rbind, lapply(1:30, function(i) {
    e <- stats::rnorm(12, 0, shock)
    x <- y <- numeric(12)
    x[1] <- stats::rnorm(1); y[1] <- stats::rnorm(1)
    for (t in 2:12) {
      x[t] <- a * x[t - 1] + e[t]
      y[t] <- a * y[t - 1] + e[t]
    }
    data.frame(id = i, time = (0:11) / 2, X = x + stats::rnorm(12, 0, sqrt(.1)),
      Y = y + stats::rnorm(12, 0, sqrt(.1)))
  }))
  m <- ctModel(type = "ct", latentNames = c("x", "y"), manifestNames = c("X", "Y"),
    LAMBDA = diag(2), CINT = matrix(0, 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    T0MEANS = matrix(0, 2, 1), silent = TRUE)
  warnings <- character()
  fit <- withCallingHandlers(suppressMessages(ctFit(dat, m, backend = "julia",
    cores = 1, verbose = 0)), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(isTRUE(fit$optim$converged))
  expect_true(any(grepl("^DIFFUSION correlations at their bound", warnings)))
  expect_false(any(grepl("Random-effect correlations", warnings)))
  diff <- fit$corbound[fit$corbound$matrix == "DIFFUSION", ]
  expect_true(all(diff$at_bound))
  expect_equal(diff$rank[1], 1)
  expect_false(any(fit$corbound$at_bound[fit$corbound$matrix == "T0VAR"]))
})
