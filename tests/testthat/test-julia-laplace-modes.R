# The inner solve's choice of mode (`_laplace_probe_modes`). Where a unit's
# integrand has more than one maximum, Newton from the origin reaches whichever
# one the origin's basin holds; probes along the found mode's soft directions,
# and from its reflection through the origin, find others, and between separate
# maxima the largest term, scored exactly as the evaluation scores it
# (`_laplace_unit_primal`), decides. A point Newton
# stopped at that climbs to a higher one without a dip is not a maximum, and
# gives way to it.
#
# The fixture is a loading whose sign the data cannot fix: y = lambda * eta
# with eta centred, so a person's likelihood depends on lambda^2 and has a
# maximum either side of zero. At a population mean of zero the origin is the
# stationary point between them -- a minimum along the loading -- which Newton
# accepts at once and the evaluation scores by the gated rule's soft direction.

.modes_sign_data <- function() {
  set.seed(3)
  n <- 12; nt <- 8
  d <- data.frame(id = rep(1:n, each = nt), time = rep(0:(nt - 1), n))
  lambda <- 1.5 * sample(c(-1, 1), n, TRUE)
  eta <- unlist(lapply(1:n, function(i) {
    e <- numeric(nt)
    e[1] <- rnorm(1)
    for (t in 2:nt) e[t] <- 0.6 * e[t - 1] + rnorm(1, 0, 0.8)
    e
  }))
  d$Y1 <- lambda[d$id] * eta + rnorm(n * nt, 0, 0.3)
  d
}

.modes_sign_model <- function() {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", id = "id",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix("lambda"),
    DRIFT = matrix(-0.5), DIFFUSION = matrix(0.8), MANIFESTVAR = matrix(0.3),
    MANIFESTMEANS = matrix(0), T0MEANS = matrix(0), T0VAR = matrix(1),
    CINT = matrix(0))))
  m$pars$indvarying <- m$pars$param %in% "lambda"
  m$pars$transform[m$pars$param %in% "lambda"] <- "param"
  m
}

test_that("the inner solve ends at a maximum, scored as the evaluation scores it", {
  skip_without_julia()
  spec <- suppressMessages(ctFit(.modes_sign_data(), .modes_sign_model(),
    backend = "julia", intoverpop = "laplace", fit = FALSE))
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  objective <- ctsem:::.ctJuliaObjective(handle, hybrid = TRUE)
  probes <- function(on) invisible(ctsem:::.ctJuliaCall(
    "ContinuousTimeSEM.ctsem_set_mode_probes!", on))
  on.exit(probes(TRUE), add = TRUE)
  evaluate <- function(x) ctJuliaEvaluate(handle, x, gradient = TRUE,
    contributions = TRUE, cores = 1L)
  # Each unit's term at the modes the last evaluation left, by the function the
  # probes compare with -- it must be the evaluation's own, the gated rule at a
  # curvature that does not factor included -- and its inner objective there
  # and at the origin.
  ctsem:::.ctJuliaEval(paste(
    "function _ctsem_test_unit_terms(o, theta)",
    "  C = ContinuousTimeSEM",
    "  Ls = C._laplace_popchols(theta, o.spec)",
    "  aws = C._laplace_workspace!(o, Float64, length(theta))",
    "  n = length(o.units.members)",
    "  g(U, u) = C._laplace_unit_objective_gradient(o, U, theta, Ls, u, aws).value",
    "  hcat([C._laplace_unit_term(o, U, theta, Ls, copy(o.modes[U]), aws) for U in 1:n],",
    "    [g(U, copy(o.modes[U])) for U in 1:n], [g(U, zeros(1)) for U in 1:n])",
    "end", sep = "\n"))
  units <- function(x) matrix(as.numeric(ctsem:::.ctJuliaCall(
    "Main._ctsem_test_unit_terms", objective, ctsem:::.ctJuliaNumericVector(x))),
    ncol = 3, dimnames = list(NULL, c("term", "g", "g_origin")))

  # Population mean zero, raw sd 0.5: every origin a stationary point. Newton
  # alone stops there; with the probes every person ends at one of their two
  # maxima, above the origin.
  x <- c(0, 0.5)
  probes(FALSE)
  off <- evaluate(x)
  u <- units(x)
  expect_equal(u[, "term"], as.numeric(off$subject_loglik), tolerance = 1e-10)
  expect_equal(u[, "g"], u[, "g_origin"], tolerance = 1e-12)
  probes(TRUE)
  ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_reset_opcounts!()")
  on <- evaluate(x)
  counts <- ctsem:::.ctJuliaGet(ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_opcounts()"))
  u <- units(x)
  expect_equal(u[, "term"], as.numeric(on$subject_loglik), tolerance = 1e-10)
  expect_true(all(u[, "g"] > u[, "g_origin"] + 0.1))
  # And the probes' counts say so: each person's better maximum came from a
  # probe's Newton solve. (At the origin the curvature need not factor, and
  # then a probe goes to Newton unscreened.)
  expect_gte(counts$probe_better, nrow(u))
  expect_gte(counts$probe_newtons, counts$probe_better)
  expect_gte(counts$probe_starts, counts$probe_newtons)

  # Away from zero the origin lies in a basin, and the probes' choice is a
  # smooth function there: the gradient is the value's.
  x <- c(0.3, 0.5)
  r <- evaluate(x)
  set.seed(1)
  v <- rnorm(2)
  v <- v / sqrt(sum(v^2))
  h <- 1e-5
  fd <- (evaluate(x + h * v)$value - evaluate(x - h * v)$value) / (2 * h)
  expect_equal(fd, sum(as.numeric(r$gradient) * v), tolerance = 1e-5)
})

# The sign fixture made wide: six loadings, each random by person, so each unit
# has six effects, and its likelihood is symmetric under flipping them all. At
# loading means of zero -- the default start's -- the origin is the stationary
# point between a person's two maxima, wherever the width. A cap on the width
# probed once left every such unit there, and the objective 4800 low.
test_that("a wide unit leaves the stationary origin for a maximum", {
  skip_without_julia()
  set.seed(3)
  n <- 10; nt <- 8; p <- 6
  d <- data.frame(id = rep(1:n, each = nt), time = rep(0:(nt - 1), n))
  signs <- sample(c(-1, 1), n, TRUE)
  eta <- unlist(lapply(1:n, function(i) {
    e <- numeric(nt)
    e[1] <- rnorm(1)
    for (t in 2:nt) e[t] <- 0.6 * e[t - 1] + rnorm(1, 0, 0.8)
    e
  }))
  for (j in 1:p) d[[paste0("Y", j)]] <- signs[d$id] * (0.6 + 0.2 * j) * eta +
    rnorm(n * nt, 0, 0.4)
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", id = "id",
    manifestNames = paste0("Y", 1:p), latentNames = "eta1",
    LAMBDA = matrix(paste0("l", 1:p), p, 1), DRIFT = matrix(-0.5),
    DIFFUSION = matrix(0.8), MANIFESTVAR = diag(0.16, p),
    MANIFESTMEANS = matrix(0, p, 1), T0MEANS = matrix(0), T0VAR = matrix(1),
    CINT = matrix(0))))
  m$pars$indvarying <- m$pars$param %in% paste0("l", 1:p)
  m$pars$transform[m$pars$param %in% paste0("l", 1:p)] <- "param"
  spec <- suppressMessages(ctFit(d, m, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  objective <- ctsem:::.ctJuliaObjective(handle, hybrid = TRUE)
  x <- numeric(ctsem:::.ctBackendNpar(spec))
  ctJuliaEvaluate(handle, x, gradient = FALSE, cores = 1L)
  ctsem:::.ctJuliaEval(paste(
    "function _ctsem_test_wide_origin(o, theta)",
    "  C = ContinuousTimeSEM",
    "  Ls = C._laplace_popchols(theta, o.spec)",
    "  aws = C._laplace_workspace!(o, Float64, length(theta))",
    "  g(U, u) = C._laplace_unit_objective_gradient(o, U, theta, Ls, u, aws).value",
    "  [g(U, copy(o.modes[U])) - g(U, zeros(length(o.modes[U])))",
    "    for U in eachindex(o.modes)]",
    "end", sep = "\n"))
  rise <- as.numeric(ctsem:::.ctJuliaCall("Main._ctsem_test_wide_origin",
    objective, ctsem:::.ctJuliaNumericVector(x)))
  expect_length(rise, n)
  expect_true(all(rise > 0.1))
})

# The loading fixture made asymmetric: the latent starts at a known mean of 1,
# so the data say the loading's sign weakly, through the first few means, and
# the sign's two maxima no longer score alike. With the population mean at 0.5
# the origin lies on the positive side, and for a person whose data favour the
# negative one Newton from the origin stops at a maximum the data pin (precision
# above `_LAPLACE_PROBE_SOFT`), where no soft direction leads anywhere.
.modes_reflect_data <- function() {
  set.seed(3)
  n <- 6; nt <- 10
  d <- data.frame(id = rep(1:n, each = nt), time = rep(0:(nt - 1), n))
  eta <- unlist(lapply(1:n, function(i) {
    e <- numeric(nt)
    e[1] <- 1 + rnorm(1)
    for (t in 2:nt) e[t] <- 0.6 * e[t - 1] + rnorm(1, 0, 0.8)
    e
  }))
  d$Y1 <- -1.2 * eta + rnorm(n * nt, 0, 0.3)
  d
}

test_that("a better mode past the origin is found where the found one is pinned", {
  skip_without_julia()
  m <- .modes_sign_model()
  m$pars$value[m$pars$matrix %in% "T0MEANS"] <- 1
  spec <- suppressMessages(ctFit(.modes_reflect_data(), m, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  objective <- ctsem:::.ctJuliaObjective(handle, hybrid = TRUE)
  x <- c(0.5, 0.5)
  r <- ctJuliaEvaluate(handle, x, gradient = FALSE, contributions = TRUE,
    cores = 1L)
  # Every maximum of each person's g_U, from a grid refined by Newton, and the
  # largest of their terms; Newton from the origin, and the precision there.
  ctsem:::.ctJuliaEval(paste(
    "function _ctsem_test_reflect_scan(o, theta)",
    "  C = ContinuousTimeSEM",
    "  Ls = C._laplace_popchols(theta, o.spec)",
    "  aws = C._laplace_workspace!(o, Float64, length(theta))",
    "  grid = collect(-6.0:0.01:6.0)",
    "  out = zeros(length(o.units.members), 4)",
    "  for U in eachindex(o.units.members)",
    "    g = [C._laplace_unit_objective_gradient(o, U, theta, Ls, [u], aws).value",
    "      for u in grid]",
    "    best = -Inf",
    "    for i in 2:length(grid)-1",
    "      (g[i] > g[i-1] && g[i] >= g[i+1]) || continue",
    "      m = C._laplace_newton_unit_mode(o, U, theta, Ls, aws, [grid[i]])",
    "      best = max(best, C._laplace_unit_primal(o, U, theta, Ls, copy(m.u), aws).term)",
    "    end",
    "    r0 = C._laplace_newton_unit_mode(o, U, theta, Ls, aws, zeros(1))",
    "    M = C._laplace_unit_curvature(o, U, theta, Ls, copy(r0.u))",
    "    out[U, :] = [best, r0.u[1], C._laplace_unit_primal(o, U, theta, Ls,",
    "      copy(r0.u), aws).term, C._laplace_block_dense(M, o.units.blocks[U], 1)[1, 1]]",
    "  end",
    "  out",
    "end", sep = "\n"))
  scan <- matrix(as.numeric(ctsem:::.ctJuliaCall("Main._ctsem_test_reflect_scan",
    objective, ctsem:::.ctJuliaNumericVector(x))), ncol = 4,
    dimnames = list(NULL, c("best", "u_origin", "term_origin", "precision")))
  # Each person's term is that of their best maximum.
  expect_equal(as.numeric(r$subject_loglik), scan[, "best"], tolerance = 1e-6)
  # And the case under test is here: a person whose origin's maximum is pinned,
  # so no soft-direction probe runs, and lower than the one past the origin.
  pinned <- scan[, "precision"] > 3.5 & scan[, "u_origin"] > 0 &
    scan[, "best"] > scan[, "term_origin"] + 0.5
  expect_true(any(pinned))
})

# Where every unit's inner problem is a concave quadratic -- Gaussian
# indicators, a linear model, effects only in mean vectors -- there is one
# maximum, and the probes are skipped (`.ctJuliaLaplaceQuadratic()`). The
# claim to hold is that nothing else changes: the same values and gradients as
# with the probes running. The sign fixture above is the other side: its
# effect is a loading, which is why it has two maxima, and it keeps them.
test_that("a model of concave units goes unprobed and evaluates as when probed", {
  skip_without_julia()
  data <- .modes_sign_data()
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct", id = "id",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix("drift"), DIFFUSION = matrix("diff"), MANIFESTVAR = matrix("mvar"),
    MANIFESTMEANS = matrix("mm"), T0MEANS = matrix("t0m"), T0VAR = matrix(1),
    CINT = matrix(0))))
  model$pars$indvarying <- model$pars$param %in% c("mm", "t0m")
  spec <- suppressMessages(ctFit(data, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  expect_true(isTRUE(spec$laplace$quadratic))
  signed <- suppressMessages(ctFit(data, .modes_sign_model(), backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  expect_false(isTRUE(signed$laplace$quadratic))

  probed <- spec
  probed$laplace$quadratic <- FALSE
  handles <- lapply(list(spec, probed), function(s)
    structure(s, class = c("ctJuliaModel", "ctFitModel")))
  # The engine holds what it was sent: a flag it ignored would pass the
  # comparison below for the wrong reason.
  field <- ctsem:::.ctJuliaEval("o -> o.quadratic")
  expect_identical(vapply(handles, function(h)
    isTRUE(field(ctsem:::.ctJuliaObjective(h, hybrid = TRUE))), logical(1)),
    c(TRUE, FALSE))
  set.seed(9)
  npar <- ctsem:::.ctBackendNpar(spec)
  for (i in 1:3) {
    x <- stats::rnorm(npar, 0, 0.5)
    a <- ctJuliaEvaluate(handles[[1]], x, gradient = TRUE, cores = 1L)
    b <- ctJuliaEvaluate(handles[[2]], x, gradient = TRUE, cores = 1L)
    expect_equal(a$value, b$value, tolerance = 1e-12)
    expect_equal(a$gradient, b$gradient, tolerance = 1e-10)
  }
})

# Between those two cases: a unit can have a second maximum only within the
# span of the loadings' rows for the effects that reach its data and are not
# the quadratic test's mean effects (`_laplace_probe_span`). Across the rest
# its likelihood is flat, or concave given the rest. The sign fixture with a
# second item, whose loading is the person's and which half the people never
# answer, and an intercept on each item with an own variance beside one shared
# dimension, on the layer (the hybrid would integrate the own variances). For
# those without the second item, its loading, intercept and own variance reach
# nothing, and the shared dimension's move of the first intercept is undone by
# that intercept's own coordinate: a flat direction that is no coordinate
# axis. The intercepts are mean effects. So the span is the loading's row,
# and nothing for those without the item. Each exclusion is compared with the
# probes without it -- the mean effects unsent, then pruning off, which counts
# every effect as reaching -- for the same values and fewer starts.
test_that("the probes search only where a unit can have a second maximum", {
  skip_without_julia()
  withr::local_options(ctsem.laplace.hybrid = FALSE)
  data <- .modes_sign_data()
  data$Y2 <- data$Y1 / 1.5 + stats::rnorm(nrow(data), 0, 0.3)
  data$Y2[data$id > 6] <- NA
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct", id = "id",
    manifestNames = c("Y1", "Y2"), latentNames = "eta1",
    LAMBDA = matrix(c("1", "lambda"), 2, 1), DRIFT = matrix(-0.5),
    DIFFUSION = matrix(0.8), MANIFESTVAR = diag(0.3, 2),
    MANIFESTMEANS = matrix(c("mm1", "mm2"), 2, 1), T0MEANS = matrix(0),
    T0VAR = matrix(1), CINT = matrix(0))))
  model$pars$indvarying <- model$pars$param %in% c("lambda", "mm1", "mm2")
  model$pars$transform[model$pars$param %in% "lambda"] <- "param"
  spec <- suppressMessages(ctFit(data, model, backend = "julia",
    intoverpop = "laplace", poprank = 1, popunique = "MANIFESTMEANS",
    fit = FALSE))
  level <- spec$laplace$levels[[1]]
  expect_null(spec$laplace$hybrid)
  expect_false(isTRUE(spec$laplace$quadratic))
  expect_equal(unname(unlist(spec$laplace$affine)), level$param != "lambda")
  full <- spec
  full$laplace$affine <- NULL
  handles <- lapply(list(spec, full), function(s)
    structure(s, class = c("ctJuliaModel", "ctFitModel")))
  # The engine holds the mean effects and the inert ones.
  held <- ctsem:::.ctJuliaEval(
    "o -> [sum(sum, o.affine; init = 0), sum(g -> sum(length, g; init = 0), o.inert; init = 0)]")
  expect_equal(as.numeric(held(ctsem:::.ctJuliaObjective(handles[[1]],
    hybrid = TRUE))), c(2, nrow(spec$laplace$inert)))

  prune <- function(on) invisible(ctsem:::.ctJuliaCall(
    "ContinuousTimeSEM.ctsem_set_laplace_prune!", on))
  on.exit(prune(TRUE), add = TRUE)
  run <- function(handle, x) {
    ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_reset_opcounts!()")
    r <- ctJuliaEvaluate(handle, x, gradient = TRUE, cores = 1L)
    list(value = r$value, gradient = r$gradient, starts = ctsem:::.ctJuliaGet(
      ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_opcounts()"))$probe_starts)
  }
  # Small population sds, so every direction is soft and probed, at two of
  # the loading's means.
  for (mean in c(0, 0.4)) {
    x <- rep(0.3, ctsem:::.ctBackendNpar(spec))
    for (f in c("sd_index", "unique_index")) x[level[[f]][level[[f]] > 0]] <- -1.5
    x[level$load_index[level$load_index > 0]] <- 0.1
    x[level$re_index[level$param == "lambda"]] <- mean
    a <- run(handles[[1]], x)
    b <- run(handles[[2]], x)
    prune(FALSE)
    c <- run(handles[[2]], x)
    prune(TRUE)
    expect_equal(a$value, c$value, tolerance = 1e-8)
    expect_equal(b$value, c$value, tolerance = 1e-8)
    expect_equal(a$gradient, c$gradient, tolerance = 1e-6)
    expect_equal(b$gradient, c$gradient, tolerance = 1e-6)
    expect_lt(a$starts, b$starts)
    expect_lt(b$starts, c$starts)
  }
})
