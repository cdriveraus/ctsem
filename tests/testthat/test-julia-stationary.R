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
.stat_model <- function(type = "ct",
    DRIFT = matrix(c("d11", "d21", "PARS[1,1]", "d22"), 2, 2),
    PARS = matrix("p1", 1, 1)) {
  # The effect is set below; ctModel warns that none is, first.
  m <- suppressWarnings(suppressMessages(ctModel(type = type, LAMBDA = diag(2),
    latentNames = c("a", "b"), manifestNames = c("y1", "y2"),
    DRIFT = DRIFT, PARS = PARS, MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(.3, 2), CINT = matrix(c("c1", "c2"), 2, 1),
    n.TIpred = 1, TIpredNames = "g", tipredDefault = FALSE, silent = TRUE)))
  m$pars$indvarying <- FALSE
  m$pars$g_effect[m$pars$param %in% "c1"] <- TRUE
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
  statedep <- .stat_model(PARS = matrix(c("p1", "p1 * b"), 2, 1),
    DRIFT = matrix(c("d11", "d21", "PARS[2,1]", "d22"), 2, 2))
  expect_error(prep(statedep), "DRIFT\\[1,2\\] depend on the latent states")

  # A latent that never settles.
  static <- .stat_model()
  rows <- static$pars$matrix %in% "DRIFT" & static$pars$row == 2
  static$pars$param[rows] <- NA
  static$pars$value[rows] <- 0
  expect_error(prep(static), "DRIFT is fixed at zero for b")

  # A single element of a latent another drives: a's DRIFT row reads PARS.
  named <- .stat_model()
  named$pars$param[named$pars$matrix %in% "T0MEANS"] <- "stationary"
  expect_error(suppressMessages(ctFit(dat, named, backend = "julia",
    fit = FALSE)), "no other latent drives, and a is driven")
})

# Single elements: a T0MEANS cell or a T0VAR diagonal cell named 'stationary'.
# Latent a is driven by nothing, b by a.
.stat_elements <- function(T0MEANS = matrix(c("stationary", "t0b"), 2, 1),
    T0VAR = matrix(c("stationary", "t0ab", 0, "t0bb"), 2, 2)) {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", LAMBDA = diag(2),
    latentNames = c("a", "b"), manifestNames = c("y1", "y2"),
    DRIFT = matrix(c("d11", "d21", 0, "d22"), 2, 2), MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(.3, 2), CINT = matrix(c("c1", "c2"), 2, 1),
    T0MEANS = T0MEANS, T0VAR = T0VAR, silent = TRUE)))
  m$pars$indvarying <- FALSE
  m
}

test_that("a 'stationary' T0 element starts only that latent's mean or variance there", {
  dat <- .stat_data()[, c("id", "time", "y1", "y2")]
  spec <- suppressMessages(ctFit(dat, .stat_elements(), backend = "julia", fit = FALSE))
  expect_false(isTRUE(spec$stationary))
  npar <- ctsem:::.ctBackendNpar(spec)
  set.seed(7)
  raw <- stats::rnorm(npar, 0, .3)
  k <- ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_kalman",
    ctsem:::.ctJuliaObjective(spec), ctsem:::.ctJuliaNumericVector(raw))
  get <- .stat_matrices(spec, raw)
  A <- get("DRIFT"); Q <- get("DIFFUSIONcov"); cint <- as.numeric(get("CINT"))
  # a starts at its own stationary mean and variance; b at its parameters; and
  # their covariance is T0VAR's correlation at a's stationary sd -- which is
  # what the summary reports too.
  expect_equal(k$eta[1, 1, 1], -cint[1] / A[1, 1], tolerance = 1e-8)
  expect_equal(k$etacov[1, 1, 1, 1], Q[1, 1] / (-2 * A[1, 1]), tolerance = 1e-8)
  expect_equal(k$eta[1, 1, 2], as.numeric(get("T0MEANS"))[2], tolerance = 1e-10)
  expect_equal(k$etacov[1, 1, , ], get("T0cov"), tolerance = 1e-10)
  expect_gt(abs(get("T0cov")[1, 2]), 1e-6)
  # The latent's own parameters are declared in PARS, which the cells read.
  pars <- spec$parameter_table
  expect_true(all(c("d11", "c1") %in% pars$param[pars$matrix %in% "PARS"]))

  # The reverse pass through the rewritten cells.
  adjoint <- as.numeric(ctJuliaEvaluate(spec, raw, gradient = TRUE,
    gradient_method = "adjoint")$gradient)
  forward <- as.numeric(ctJuliaEvaluate(spec, raw, gradient = TRUE,
    gradient_method = "forward")$gradient)
  expect_equal(adjoint, forward, tolerance = 1e-9)

  # Only what is named: the mean alone leaves a's variance a parameter.
  meanonly <- suppressMessages(ctFit(dat, .stat_elements(T0VAR = matrix(c("t0aa", "t0ab", 0,
    "t0bb"), 2, 2)), backend = "julia", fit = FALSE))
  names <- as.character(ctsem:::.ctBackendRawParameterNames(list(model_spec = meanonly),
    ctsem:::.ctBackendNpar(meanonly)))
  expect_true("t0aa" %in% names)
  expect_false(any(c("T0mean_a", "t0a") %in% names))
})

test_that("a single latent with both elements named starts as stationary = TRUE does", {
  dat <- .stat_data()[, c("id", "time", "y1")]
  one <- function(T0MEANS, T0VAR) {
    m <- suppressWarnings(suppressMessages(ctModel(type = "ct", LAMBDA = diag(1),
      latentNames = "a", manifestNames = "y1", DRIFT = matrix("d11"),
      DIFFUSION = matrix("q11"), CINT = matrix("c1"), MANIFESTMEANS = matrix(0),
      MANIFESTVAR = matrix(.3), T0MEANS = T0MEANS, T0VAR = T0VAR, silent = TRUE)))
    m$pars$indvarying <- FALSE
    m
  }
  elements <- suppressMessages(ctFit(dat, one(matrix("stationary"), matrix("stationary")),
    backend = "julia", fit = FALSE))
  whole <- suppressMessages(ctFit(dat, one(matrix("t0m"), matrix("t0v")), backend = "julia",
    fit = FALSE, stationary = TRUE))
  named <- function(spec) as.character(ctsem:::.ctBackendRawParameterNames(
    list(model_spec = spec), ctsem:::.ctBackendNpar(spec)))
  ne <- named(elements); nw <- named(whole)
  expect_setequal(ne, nw)
  set.seed(8)
  raw <- setNames(stats::rnorm(length(nw), 0, .3), nw)
  expect_equal(ctJuliaEvaluate(elements, raw[ne], gradient = FALSE)$value,
    ctJuliaEvaluate(whole, raw[nw], gradient = FALSE)$value, tolerance = 1e-10)
})

test_that("'stationary' elements are refused where they are not single latents' own", {
  dat <- .stat_data(4)[, c("id", "time", "y1", "y2")]
  prep <- function(m, ...) suppressMessages(ctFit(dat, m, backend = "julia", fit = FALSE, ...))
  # b is driven by a.
  expect_error(prep(.stat_elements(T0MEANS = matrix(c("t0a", "stationary"), 2, 1))),
    "b is driven through DRIFT\\[2,1\\]")
  # An off-diagonal covariance.
  expect_error(prep(.stat_elements(T0VAR = matrix(c("t0aa", "stationary", 0, "t0bb"), 2, 2))),
    "Only T0MEANS cells and the diagonal of T0VAR")
  expect_error(suppressMessages(ctFit(dat, .stat_elements(), backend = "stan", fit = FALSE)),
    "needs backend = 'julia'")
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
