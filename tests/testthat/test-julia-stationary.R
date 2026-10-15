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

  # A single stationary element in a latent that depends on one that never
  # settles: a is driven by b, a random walk.
  walk <- .stat_model()
  rows <- walk$pars$matrix %in% "DRIFT" & walk$pars$row == 2
  walk$pars$param[rows] <- NA
  walk$pars$value[rows] <- 0
  walk$pars$param[walk$pars$matrix %in% "T0MEANS" & walk$pars$row == 1] <- "stationary"
  expect_error(suppressMessages(ctFit(dat, walk, backend = "julia", fit = FALSE)),
    "DRIFT is fixed at zero for b")
})

# Single elements: T0MEANS and T0VAR cells named 'stationary', on the model
# above -- a driven by b through PARS, and a TI predictor on a's CINT, so the
# stationary moments are the joint system's and differ by subject.
.stat_marked <- function(T0MEANS = c("stationary", "t0b"),
    T0VAR = matrix(c("stationary", "t0ab", 0, "t0bb"), 2, 2), type = "ct") {
  m <- .stat_model(type = type)
  for (cell in list(list("T0MEANS", matrix(T0MEANS, 2, 1)), list("T0VAR", T0VAR))) {
    for (r in 1:2) for (c in 1:ncol(cell[[2]])) {
      k <- m$pars$matrix %in% cell[[1]] & m$pars$row == r & m$pars$col == c
      if (!any(k)) next
      value <- cell[[2]][r, c]
      if (value %in% "0") {
        m$pars$param[k] <- NA
        m$pars$value[k] <- 0
      } else {
        m$pars$param[k] <- value
        m$pars$value[k] <- NA
      }
    }
  }
  m
}

.stat_names <- function(spec) as.character(ctsem:::.ctBackendRawParameterNames(
  list(model_spec = spec), ctsem:::.ctBackendNpar(spec)))

test_that("a 'stationary' element takes the joint system's stationary value, and only it", {
  dat <- .stat_data()
  spec <- suppressMessages(ctFit(dat, .stat_marked(), backend = "julia", fit = FALSE))
  expect_false(isTRUE(spec$stationary))
  free <- suppressMessages(ctFit(dat, .stat_marked(T0MEANS = c("t0a", "t0b"),
    T0VAR = matrix(c("t0aa", "t0ab", 0, "t0bb"), 2, 2)), backend = "julia", fit = FALSE))
  nm <- .stat_names(spec)
  nf <- .stat_names(free)
  # The marked cells are no longer parameters; nothing else changed.
  expect_setequal(setdiff(nf, nm), c("t0a", "t0aa"))
  set.seed(7)
  rawf <- setNames(stats::rnorm(length(nf), 0, .3), nf)
  raw <- unname(rawf[nm])
  rawf <- unname(rawf)

  objective <- ctsem:::.ctJuliaObjective(spec)
  k <- ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_kalman", objective,
    ctsem:::.ctJuliaNumericVector(raw))
  kf <- ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_kalman",
    ctsem:::.ctJuliaObjective(free), ctsem:::.ctJuliaNumericVector(rawf))
  first <- match(unique(k$subject), k$subject)
  for (i in seq_along(first)) {
    get <- .stat_matrices(spec, raw, tipreds = spec$tipred_data[i, ])
    expected <- .stat_moments(get("DRIFT"), get("CINT"), get("DIFFUSIONcov"))
    f <- first[i]
    # a: the joint system's stationary mean and variance.
    expect_equal(k$eta[1, f, 1], expected$mean[1], tolerance = 1e-8)
    expect_equal(k$etacov[1, f, 1, 1], expected$cov[1, 1], tolerance = 1e-8)
    # b as written, and T0VAR's correlation between them.
    expect_equal(k$eta[1, f, 2], kf$eta[1, f, 2], tolerance = 1e-10)
    expect_equal(k$etacov[1, f, 2, 2], kf$etacov[1, f, 2, 2], tolerance = 1e-10)
    expect_equal(k$etacov[1, f, 1, 2] / sqrt(k$etacov[1, f, 1, 1] * k$etacov[1, f, 2, 2]),
      kf$etacov[1, f, 1, 2] / sqrt(kf$etacov[1, f, 1, 1] * kf$etacov[1, f, 2, 2]),
      tolerance = 1e-4)
    # The summary reports the start the filter used.
    expect_equal(as.numeric(get("T0MEANS")), k$eta[1, f, ], tolerance = 1e-8)
    expect_equal(get("T0cov"), k$etacov[1, f, , ], tolerance = 1e-8)
  }
  expect_gt(abs(k$etacov[1, 1, 1, 2]), 1e-6)

  adjoint <- as.numeric(ctJuliaEvaluate(spec, raw, gradient = TRUE,
    gradient_method = "adjoint")$gradient)
  forward <- as.numeric(ctJuliaEvaluate(spec, raw, gradient = TRUE,
    gradient_method = "forward")$gradient)
  expect_equal(adjoint, forward, tolerance = 1e-9)
})

test_that("every element named is stationary = TRUE, and a named block is the stationary block", {
  dat <- .stat_data()
  all <- .stat_marked(T0MEANS = c("stationary", "stationary"),
    T0VAR = matrix(c("stationary", "stationary", 0, "stationary"), 2, 2))
  spec <- suppressMessages(ctFit(dat, all, backend = "julia", fit = FALSE))
  whole <- suppressMessages(ctFit(dat, .stat_model(), backend = "julia", fit = FALSE,
    stationary = TRUE))
  ns <- .stat_names(spec)
  nw <- .stat_names(whole)
  expect_setequal(ns, nw)
  set.seed(9)
  raw <- setNames(stats::rnorm(length(nw), 0, .3), nw)
  rs <- unname(raw[ns])
  expect_equal(ctJuliaEvaluate(spec, rs, gradient = FALSE)$value,
    ctJuliaEvaluate(whole, unname(raw), gradient = FALSE)$value, tolerance = 1e-9)

  adjoint <- as.numeric(ctJuliaEvaluate(spec, rs, gradient = TRUE,
    gradient_method = "adjoint")$gradient)
  forward <- as.numeric(ctJuliaEvaluate(spec, rs, gradient = TRUE,
    gradient_method = "forward")$gradient)
  expect_equal(adjoint, forward, tolerance = 1e-9)
  objective <- ctsem:::.ctJuliaObjective(spec)
  over_reverse <- ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_hessian",
    objective, ctsem:::.ctJuliaNumericVector(rs))
  over_forward <- ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_hessian_forward",
    objective, ctsem:::.ctJuliaNumericVector(rs))
  expect_equal(as.numeric(over_reverse), as.numeric(over_forward), tolerance = 1e-8)
})

test_that("a 'stationary' element needs only its own latents to settle", {
  # b a random walk that a does not depend on: stationary = TRUE has no
  # distribution to start from, while a's elements do.
  dat <- .stat_data()
  m <- .stat_marked(T0MEANS = c("stationary", "t0b"))
  m$pars$param[m$pars$matrix %in% "DRIFT" & m$pars$row == 1 & m$pars$col == 2] <- NA
  m$pars$value[m$pars$matrix %in% "DRIFT" & m$pars$row == 1 & m$pars$col == 2] <- 0
  rows <- m$pars$matrix %in% "DRIFT" & m$pars$row == 2
  m$pars$param[rows] <- NA
  m$pars$value[rows] <- 0
  m$pars$g_effect[] <- FALSE
  spec <- suppressMessages(ctFit(dat, m, backend = "julia", fit = FALSE))
  set.seed(10)
  raw <- stats::rnorm(ctsem:::.ctBackendNpar(spec), 0, .3)
  k <- ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_kalman",
    ctsem:::.ctJuliaObjective(spec), ctsem:::.ctJuliaNumericVector(raw))
  get <- .stat_matrices(spec, raw)
  A <- get("DRIFT"); Q <- get("DIFFUSIONcov"); cint <- as.numeric(get("CINT"))
  expect_equal(k$eta[1, 1, 1], -cint[1] / A[1, 1], tolerance = 1e-8)
  expect_equal(k$etacov[1, 1, 1, 1], Q[1, 1] / (-2 * A[1, 1]), tolerance = 1e-8)
  expect_error(suppressMessages(ctFit(dat, .stat_model(DRIFT = matrix(c("d11", 0, 0, 0),
    2, 2)), backend = "julia", fit = FALSE, stationary = TRUE)), "DRIFT is fixed at zero")
})

test_that("'stationary' elements are refused where they name no stationary value", {
  dat <- .stat_data(4)
  prep <- function(m, ...) suppressMessages(ctFit(dat, m, backend = "julia", fit = FALSE, ...))
  # A covariance without both variances.
  expect_error(prep(.stat_marked(T0VAR = matrix(c("t0aa", "stationary", 0, "t0bb"), 2, 2))),
    "needs both their variances")
  # Not an initial condition.
  drift <- .stat_marked()
  drift$pars$param[drift$pars$matrix %in% "DRIFT" & drift$pars$row == 2 &
    drift$pars$col == 1] <- "stationary"
  expect_error(prep(drift), "Only T0MEANS and T0VAR cells")
  # A random initial mean with a stationary variance.
  random <- .stat_marked(T0MEANS = c("t0a", "t0b"))
  random$pars$indvarying[random$pars$param %in% "t0a"] <- TRUE
  expect_error(prep(random), "Name the T0MEANS cell 'stationary' too")
  expect_error(suppressMessages(ctFit(dat, .stat_marked(), backend = "stan",
    fit = FALSE)), "needs backend = 'julia'")
  expect_error(prep(.stat_marked(type = "dt")), "continuous time")
})

test_that("generation starts 'stationary' elements where a fit does", {
  # Fixed: the first observation is the stationary distribution, mean
  # -CINT / DRIFT = 4 and variance DIFFUSION^2 / (-2 DRIFT) = 1, plus the
  # measurement error's 0.01.
  one <- suppressWarnings(suppressMessages(ctModel(type = "ct", LAMBDA = diag(1),
    latentNames = "a", manifestNames = "y1", DRIFT = matrix(-0.5),
    DIFFUSION = matrix(1), CINT = matrix(2), MANIFESTMEANS = matrix(0),
    MANIFESTVAR = matrix(.1), T0MEANS = matrix("stationary"),
    T0VAR = matrix("stationary"), silent = TRUE)))
  set.seed(11)
  dat <- suppressMessages(ctGenerate(one, n = 4000, Tpoints = 2))
  first <- dat[dat[, "time"] == 0, "y1"]
  expect_lt(abs(mean(first) - 4), 0.1)
  expect_lt(abs(stats::var(first) - 1.01), 0.1)
  expect_error(ctGenerate(one, n = 5, Tpoints = 2, backend = "r"),
    "cannot place")
})

test_that("random effects on the dynamics take the route that keeps them fixed per subject", {
  dat <- .stat_data(6)
  m <- .stat_model()
  m$pars$indvarying[m$pars$param %in% "c1"] <- TRUE
  expect_error(suppressMessages(ctFit(dat, m, backend = "julia", fit = FALSE,
    stationary = TRUE, intoverpop = "augmented")),
    "carries the individual variation in CINT\\[1,1\\] as a state")
  expect_message(spec <- ctFit(dat, m, backend = "julia", fit = FALSE,
    stationary = TRUE), "chose 'laplace'.*stationary initial conditions")
  expect_false(is.null(spec$laplace))
  expect_true(spec$stationary)
  # Without stationarity the same model stays on the augmented route.
  plain <- suppressMessages(ctFit(dat, m, backend = "julia", fit = FALSE))
  expect_true(is.null(plain$laplace))
})
