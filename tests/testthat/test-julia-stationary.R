# ctFit(stationary = TRUE): each subject's latent processes start from the
# distribution DRIFT, CINT and DIFFUSION settle into, in place of T0MEANS and
# T0VAR. The moments are checked against that distribution worked out here by
# a Kronecker solve -- not by the engine's Lyapunov solver or its intercept
# solve -- at the filter's first prior, per subject when a TI predictor moves
# CINT, and in what the summary reports. The reverse pass is checked against
# forward mode, and each refusal by the setting or cell that causes it.

skip_without_julia()

# Two latents, a DRIFT cell written by a PARS transform (so the predict group
# has to run before the stationary moments read it), and a TI predictor on
# CINT (so subjects have different stationary means).
.stat_model <- function(...) {
  m <- suppressMessages(ctModel(type = "ct", LAMBDA = diag(2),
    latentNames = c("a", "b"), manifestNames = c("y1", "y2"),
    DRIFT = matrix(c("d11", "d21", "PARS[1,1]", "d22"), 2, 2),
    PARS = matrix("p1", 1, 1), MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(.3, 2), CINT = matrix(c("c1", "c2"), 2, 1),
    n.TIpred = 1, TIpredNames = "g", tipredDefault = FALSE, silent = TRUE,
    ...))
  m$pars$indvarying <- FALSE
  m
}

.stat_data <- function(nsubjects = 8, seed = 3) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(nsubjects), function(i) {
    time <- cumsum(c(0, stats::rexp(3, 1)))
    data.frame(id = i, time = time, y1 = stats::rnorm(4), y2 = stats::rnorm(4),
      g = i / nsubjects - 0.5)
  }))
}

# The stationary moments by Kronecker solve.
.stat_moments <- function(A, cint, Q) {
  n <- nrow(A)
  list(mean = -solve(A, cint), cov = matrix(-solve(kronecker(diag(n), A) +
    kronecker(A, diag(n)), as.numeric(Q)), n, n))
}

# One materialised matrix from the engine's summary layout.
.stat_matrices <- function(spec, raw, tipreds = NULL) {
  objective <- ctsem:::.ctJuliaObjective(spec)
  layout <- ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_parameter_layout",
    objective)
  arguments <- list("ContinuousTimeSEM.ctsem_parameter_matrices", objective,
    ctsem:::.ctJuliaNumericVector(raw))
  if (length(tipreds)) arguments$tipreds <- ctsem:::.ctJuliaVector(tipreds)
  flat <- as.numeric(do.call(ctsem:::.ctJuliaCall, arguments))
  function(name) {
    i <- which(unlist(layout$matrix) == name)
    nr <- unlist(layout$nrow)[i]
    matrix(flat[unlist(layout$offset)[i] + seq_len(nr * unlist(layout$ncol)[i])],
      nr)
  }
}

test_that("each subject starts from the stationary distribution of its own dynamics", {
  m <- .stat_model()
  m$pars$g_effect[m$pars$param %in% "c1"] <- TRUE
  dat <- .stat_data()
  spec <- suppressMessages(ctFit(dat, m, backend = "julia", fit = FALSE,
    stationary = TRUE))
  expect_true(spec$stationary)
  t0 <- spec$parameter_table$matrix %in% c("T0MEANS", "T0VAR")
  expect_true(all(is.na(spec$parameter_table$parnumber[t0])))

  npar <- ctsem:::.ctBackendNpar(spec)
  set.seed(5)
  raw <- stats::rnorm(npar, 0, .3)
  k <- ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_kalman",
    ctsem:::.ctJuliaObjective(spec), ctsem:::.ctJuliaNumericVector(raw))
  first <- match(unique(k$subject), k$subject)
  means <- vector("list", length(first))
  for (i in seq_along(first)) {
    get <- .stat_matrices(spec, raw, tipreds = spec$tipred_data[i, ])
    expected <- .stat_moments(get("DRIFT"), get("CINT"), get("DIFFUSIONcov"))
    expect_equal(k$eta[1, first[i], ], as.numeric(expected$mean),
      tolerance = 1e-8)
    expect_equal(k$etacov[1, first[i], , ], expected$cov, tolerance = 1e-8)
    # The summary reports the same start the filter used.
    expect_equal(as.numeric(get("T0MEANS")), as.numeric(expected$mean),
      tolerance = 1e-8)
    expect_equal(get("T0cov"), expected$cov, tolerance = 1e-8)
    means[[i]] <- expected$mean
  }
  # The TI predictor really did move the start, or the loop above proved
  # nothing per subject.
  expect_gt(abs(means[[1]][1] - means[[length(means)]][1]), 1e-3)
})

test_that("the reverse pass differentiates the stationary start", {
  m <- .stat_model()
  m$pars$g_effect[m$pars$param %in% "c1"] <- TRUE
  spec <- suppressMessages(ctFit(.stat_data(), m, backend = "julia",
    fit = FALSE, stationary = TRUE))
  npar <- ctsem:::.ctBackendNpar(spec)
  set.seed(6)
  raw <- stats::rnorm(npar, 0, .3)
  adjoint <- as.numeric(ctJuliaEvaluate(spec, raw, gradient = TRUE,
    gradient_method = "adjoint")$gradient)
  forward <- as.numeric(ctJuliaEvaluate(spec, raw, gradient = TRUE,
    gradient_method = "forward")$gradient)
  expect_equal(adjoint, forward, tolerance = 1e-9)
  # Every parameter the start reads has a gradient through it, so a dropped
  # term would show; and none is exactly zero by accident of the data.
  expect_true(all(abs(adjoint) > 0))
  objective <- ctsem:::.ctJuliaObjective(spec)
  over_reverse <- ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_hessian",
    objective, ctsem:::.ctJuliaNumericVector(raw))
  over_forward <- ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_hessian_forward",
    objective, ctsem:::.ctJuliaNumericVector(raw))
  expect_equal(as.numeric(over_reverse), as.numeric(over_forward),
    tolerance = 1e-8)
})

test_that("stationary = TRUE is refused where there is no one distribution to start from", {
  dat <- .stat_data(4)
  prep <- function(m, ...) suppressMessages(ctFit(dat, m, backend = "julia",
    fit = FALSE, stationary = TRUE, ...))
  expect_error(suppressMessages(ctFit(dat, .stat_model(), backend = "stan",
    fit = FALSE, stationary = TRUE)), "needs backend = 'julia'")
  expect_error(prep(.stat_model(type = "dt")), "continuous time")
  expect_error(suppressMessages(ctFit(dat, .stat_model(), backend = "julia",
    fit = FALSE, stationary = NA)), "TRUE or FALSE")

  # A state reaching DRIFT through PARS.
  statedep <- .stat_model()
  statedep$pars$param[statedep$pars$matrix %in% "PARS"] <- "p1 * b"
  expect_error(prep(statedep), "DRIFT\\[1,2\\] depend on the latent states")

  # A latent that never settles.
  static <- .stat_model()
  rows <- static$pars$matrix %in% "DRIFT" & static$pars$row == 2
  static$pars$param[rows] <- NA
  static$pars$value[rows] <- 0
  expect_error(prep(static), "DRIFT is fixed at zero for b")

  # The per-cell spelling the option once had.
  named <- .stat_model()
  named$pars$param[named$pars$matrix %in% "T0MEANS"] <- "stationary"
  expect_error(suppressMessages(ctFit(dat, named, backend = "julia",
    fit = FALSE)), "no longer supported")
})

test_that("random effects on the dynamics take the route that keeps them fixed per subject", {
  dat <- .stat_data(6)
  m <- .stat_model()
  m$pars$indvarying[m$pars$param %in% "c1"] <- TRUE
  expect_error(suppressMessages(ctFit(dat, m, backend = "julia", fit = FALSE,
    stationary = TRUE, intoverpop = "augmented")),
    "carries the individual variation in CINT\\[1,1\\] as a state")
  expect_message(spec <- ctFit(dat, m, backend = "julia", fit = FALSE,
    stationary = TRUE), "chose 'laplace'.*stationary=TRUE")
  expect_false(is.null(spec$laplace))
  expect_true(spec$stationary)
  # Without stationarity the same model stays on the augmented route.
  plain <- suppressMessages(ctFit(dat, m, backend = "julia", fit = FALSE))
  expect_true(is.null(plain$laplace))
})
