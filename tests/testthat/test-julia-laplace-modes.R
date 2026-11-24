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
# accepts at once and the evaluation scores by the soft-direction quadrature.

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
  objective <- ctsem:::.ctJuliaObjective(handle)
  probes <- function(on) invisible(ctsem:::.ctJuliaCall(
    "ContinuousTimeSEM.ctsem_set_mode_probes!", on))
  on.exit(probes(TRUE), add = TRUE)
  # Which maximum a unit ends at, scored alone: the combination of a unit's
  # maxima (`_laplace_mode_sum`) is the next test's.
  set_sum <- function(on) invisible(ctsem:::.ctJuliaCall(
    "ContinuousTimeSEM.ctsem_set_mode_sum!", on))
  set_sum(FALSE)
  on.exit(set_sum(TRUE), add = TRUE)
  evaluate <- function(x) ctJuliaEvaluate(handle, x, gradient = TRUE,
    contributions = TRUE, cores = 1L)
  # Each unit's term at the modes the last evaluation left, by the function the
  # probes compare with -- it must be the evaluation's own, the soft-direction quadrature at a
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
  # smooth function there, as is the combination of each unit's maxima: the gradient
  # is the value's.
  set_sum(TRUE)
  x <- c(0.3, 0.5)
  r <- evaluate(x)
  set.seed(1)
  v <- rnorm(2)
  v <- v / sqrt(sum(v^2))
  h <- 1e-5
  fd <- (evaluate(x + h * v)$value - evaluate(x - h * v)$value) / (2 * h)
  expect_equal(fd, sum(as.numeric(r$gradient) * v), tolerance = 1e-5)
})

# A unit with maxima a valley parts is scored by the best one's term, handed to
# the next as it fades (`_laplace_mode_sum`, `_laplace_mode_shares`), which
# adds nothing a single maximum would not: at a population mean of zero each
# person's two maxima on the sign fixture are mirror images, and the
# combination is either's term. Away from zero one is the better, and where it
# is at full weight its term is the unit's whatever the search found beside it:
# averaged by their masses, a lower maximum took a share, and the term moved
# wherever a search found or missed one. A mode's smallest eigenvalue, which
# fades it out at a fold, has the gradient its differences give with the mode
# solved again (`_laplace_mode_fade`); the combined objective's own gradient
# against its differences is the test above.
test_that("a unit with two maxima is scored by the better one's term", {
  skip_without_julia()
  spec <- suppressMessages(ctFit(.modes_sign_data(), .modes_sign_model(),
    backend = "julia", intoverpop = "laplace", fit = FALSE))
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  objective <- ctsem:::.ctJuliaObjective(handle)
  set_sum <- function(on) invisible(ctsem:::.ctJuliaCall(
    "ContinuousTimeSEM.ctsem_set_mode_sum!", on))
  on.exit(set_sum(TRUE), add = TRUE)
  evaluate <- function(x) as.numeric(ctJuliaEvaluate(handle, x, gradient = FALSE,
    contributions = TRUE, cores = 1L)$subject_loglik)
  x <- c(0, 0.5)
  set_sum(FALSE)
  one <- evaluate(x)
  set_sum(TRUE)
  both <- evaluate(x)
  expect_equal(both, one, tolerance = 1e-8)
  x <- c(0.3, 0.5)
  set_sum(FALSE)
  one <- evaluate(x)
  set_sum(TRUE)
  both <- evaluate(x)
  ctsem:::.ctJuliaEval(
    "_ctsem_test_others(o) = Float64[length(a) for a in o.alt_modes]")
  others <- as.numeric(ctsem:::.ctJuliaCall("Main._ctsem_test_others", objective))
  expect_gte(sum(others > 0), 3)
  expect_equal(both, one, tolerance = 1e-10)
  ctsem:::.ctJuliaEval(paste(
    "function _ctsem_test_fade(o, theta)",
    "  C = ContinuousTimeSEM",
    "  C.ctsem_laplace_evaluate(o, theta; gradient = false)",
    "  U = 1",
    "  u = copy(o.modes[U])",
    "  f = C._laplace_mode_fade(o, U, theta, u; gradient = true)",
    "  aws = C._laplace_workspace!(o, Float64, length(theta))",
    "  lam(x) = C._laplace_mode_fade(o, U, x, copy(C._laplace_newton_unit_mode(o, U, x,",
    "    C._laplace_popchols(x, o.spec), aws, u).u)).lam",
    "  h = 1e-5",
    "  e(j) = [k == j ? h : 0.0 for k in eachindex(theta)]",
    "  fd = [(lam(theta .+ e(j)) - lam(theta .- e(j))) / (2h) for j in eachindex(theta)]",
    "  vcat(f.lam, f.grad, fd)",
    "end", sep = "\n"))
  g <- as.numeric(ctsem:::.ctJuliaCall("Main._ctsem_test_fade", objective,
    ctsem:::.ctJuliaNumericVector(c(0.3, 0.5))))
  expect_gt(g[1], 0)
  expect_equal(g[2:3], g[4:5], tolerance = 1e-5)
})

# Three binary items on one latent, a random intercept and a random auto-effect
# under -log1p_exp: the bench's gD1 (dev/optimbench/cells.R), 60 people at 8
# occasions. At the estimate below 11 people have two maxima under one
# curvature sd apart, with a valley of hundredths of a nat between them, where
# the valley's fade is part way. Summed as masses, those maxima added up to
# 0.5 nats to each one's term, and a fit grew more of them by widening the
# random drift's sd, 3.5 exact nats down. Where the best maximum's fade is part
# way the next takes the rest of its share. The combination is for optimising,
# not an exact integral: the test is that it never scores a unit above its
# best single term, and that its gradient, with the valleys' fades moving, is
# its value's.
.modes_binary_data <- function() {
  set.seed(1)
  rows <- lapply(seq_len(60), function(i) {
    raw <- stats::rnorm(1, 1, 1.5); a <- -log1p(exp(-raw)); decay <- exp(a)
    baseline <- stats::rnorm(1, 0, 1)
    intercept <- (baseline / a) * (decay - 1)
    innovation <- sqrt((exp(2 * a) - 1) / (2 * a))
    latent <- numeric(8); latent[1] <- stats::rnorm(1, 0, 0.7)
    for (t in 1:7) latent[t + 1L] <- decay * latent[t] + intercept +
      stats::rnorm(1, 0, innovation)
    Y <- sapply(c(0, 0.5, -0.5), function(m) stats::rbinom(8, 1, stats::plogis(latent + m)))
    data.frame(id = i, time = 0:7, B1 = Y[, 1], B2 = Y[, 2], B3 = Y[, 3])
  })
  do.call(rbind, rows)
}

.modes_binary_model <- function() {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = c("B1", "B2", "B3"), latentNames = "eta1",
    LAMBDA = matrix(1, 3, 1), manifesttype = c(1, 1, 1),
    MANIFESTMEANS = matrix(c(0, "thr2", "thr3"), 3, 1), MANIFESTVAR = matrix(0, 3, 3),
    CINT = matrix("cint||TRUE"), DRIFT = matrix("drift|-log1p_exp(-param)|TRUE"),
    T0VAR = matrix(0.5), silent = TRUE)))
  m$pars$indvarying <- m$pars$param %in% c("drift", "cint")
  m
}

test_that("maxima that share a basin add nothing to a unit's term", {
  skip_without_julia()
  spec <- suppressMessages(ctFit(.modes_binary_data(), .modes_binary_model(),
    backend = "julia", intoverpop = "laplace", fit = FALSE))
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  set_sum <- function(on) invisible(ctsem:::.ctJuliaCall(
    "ContinuousTimeSEM.ctsem_set_mode_sum!", on))
  on.exit(set_sum(TRUE), add = TRUE)
  # gD1's s6 fit's estimate on the bench, the combination off (localcurv-e0aae396).
  x <- c(0.0301148912, 4.5210688928, -1.1069526929, 0.0234851509, -0.0632376482,
    -0.0056830083, 2.0234889501, -1.0443493528, -0.0565364719)
  evaluate <- function(x, gradient = FALSE) ctJuliaEvaluate(handle, x,
    gradient = gradient, contributions = TRUE, cores = 1L)
  set_sum(FALSE)
  one <- as.numeric(evaluate(x)$subject_loglik)
  set_sum(TRUE)
  r <- evaluate(x, gradient = TRUE)
  both <- as.numeric(r$subject_loglik)
  changed <- which(abs(both - one) > 1e-6)
  expect_gte(length(changed), 3)
  # Without the combination a unit scores its larger term.
  expect_true(all(both <= one + 1e-10))
  # The shares and the valleys' fades all enter the gradient here, where each
  # share is inside (0, 1).
  set.seed(2)
  v <- rnorm(length(x))
  h <- 1e-5
  fd <- (evaluate(x + h * v)$value - evaluate(x - h * v)$value) / (2 * h)
  expect_equal(fd, sum(as.numeric(r$gradient) * v), tolerance = 1e-5)
})

# An anchored run follows the maxima it started with, so one born on its way
# is the next search's to find; a run takes it where it stops and goes on,
# rather than a new run starting (`ctsem_optimize`).
test_that("a search where an anchored run stops joins the maxima it missed", {
  skip_without_julia()
  spec <- suppressMessages(ctFit(.modes_binary_data(), .modes_binary_model(),
    backend = "julia", intoverpop = "laplace", fit = FALSE))
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  objective <- ctsem:::.ctJuliaObjective(handle)
  x <- c(0.0301148912, 4.5210688928, -1.1069526929, 0.0234851509, -0.0632376482,
    -0.0056830083, 2.0234889501, -1.0443493528, -0.0565364719)
  ctsem:::.ctJuliaEval(paste(
    "function _ctsem_test_take(o, x)",
    "  C = ContinuousTimeSEM",
    "  value(o) = C.ctsem_laplace_evaluate(o, x; gradient = false).value",
    "  C._ctsem_anchor_start!(o, x)",
    "  full = value(o)",
    "  several = count(!isempty, o.alt_anchors)",
    "  known0 = sum(length, o.known)",
    "  o.alt_anchors = [Vector{Float64}[] for _ in o.alt_anchors]",
    "  held = value(o)",
    "  found = C.CTSEMFoundModes()",
    "  t = C._ctsem_take_search_modes!(o, x, found)",
    "  after = value(o)",
    "  again = C._ctsem_take_search_modes!(o, x, found)",
    "  jumped = C._ctsem_take_jump_modes!(o, x, (units = found.units[1:1],",
    "    point_units = found.units[1:1], point_values = reshape(found.modes[1], :, 1)), found)",
    "  C.ctsem_laplace_clear_anchors!(o)",
    "  [full, held, t.change, t.units, after, again.units, several,",
    "    length(found.units), sum(length, o.known) - known0, jumped]",
    "end", sep = "\n"))
  r <- as.numeric(ctsem:::.ctJuliaCall("Main._ctsem_test_take", objective,
    ctsem:::.ctJuliaNumericVector(x)))
  names(r) <- c("full", "held", "change", "units", "after", "again", "several",
    "found", "known", "jumped")
  # Held to one maximum each -- the one of largest weight, which is not always
  # of largest term -- the units with several score otherwise.
  expect_gte(r[["several"]], 3)
  expect_gt(abs(r[["held"]] - r[["full"]]), 1e-3)
  # The search gives each its maxima back, and says what that moved.
  expect_equal(r[["units"]], r[["several"]])
  expect_equal(r[["change"]], r[["full"]] - r[["held"]], tolerance = 1e-8)
  expect_equal(r[["after"]], r[["full"]], tolerance = 1e-8)
  # What it joined is known from then on, once; a second search adds nothing.
  expect_gte(r[["found"]], r[["several"]])
  expect_equal(r[["known"]], r[["found"]])
  expect_equal(r[["again"]], 0)
  # A jump's mode already known adds nothing for the run to go on with.
  expect_equal(r[["jumped"]], 0)
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
  objective <- ctsem:::.ctJuliaObjective(handle)
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
  objective <- ctsem:::.ctJuliaObjective(handle)
  x <- c(0.5, 0.5)
  # Scored alone, which maximum the search reached (the combination of them
  # is `_laplace_mode_sum`'s, tested below).
  invisible(ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_mode_sum!", FALSE))
  r <- tryCatch(ctJuliaEvaluate(handle, x, gradient = FALSE, contributions = TRUE,
    cores = 1L), finally = invisible(ctsem:::.ctJuliaCall(
      "ContinuousTimeSEM.ctsem_set_mode_sum!", TRUE)))
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
# effect is a loading, which is why it has two maxima, and it keeps them. On
# the layer, where the effects are coordinates; the hybrid would integrate
# them, and neither side would probe anything. Run, the probes find every line
# the found mode's quadratic by its value, and spend no gradient.
test_that("a model of concave units goes unprobed and evaluates as when probed", {
  skip_without_julia()
  withr::local_options(ctsem.laplace.hybrid = FALSE)
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
    isTRUE(field(ctsem:::.ctJuliaObjective(h))), logical(1)),
    c(TRUE, FALSE))
  set.seed(9)
  npar <- ctsem:::.ctBackendNpar(spec)
  for (i in 1:3) {
    x <- stats::rnorm(npar, 0, 0.5)
    a <- ctJuliaEvaluate(handles[[1]], x, gradient = TRUE, cores = 1L)
    ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_reset_opcounts!()")
    b <- ctJuliaEvaluate(handles[[2]], x, gradient = TRUE, cores = 1L)
    counts <- ctsem:::.ctJuliaGet(ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_opcounts()"))
    expect_equal(a$value, b$value, tolerance = 1e-12)
    expect_equal(a$gradient, b$gradient, tolerance = 1e-10)
    expect_gt(counts$probe_values, 0)
    expect_equal(counts$probe_gradients, 0)
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
  expect_equal(as.numeric(held(ctsem:::.ctJuliaObjective(handles[[1]]))), c(2, nrow(spec$laplace$inert)))

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

# Persons nested in studies, as in a meta-analysis of panels: person-level
# random intercepts and initial means, a study-level random drift. On the layer
# a study's unit holds every person's effects, and its mode search -- Newton's
# kept curvature, the probes a block at a time (`_laplace_probe_nested`), the
# partner modes' fold test -- works from the unit's blocks
# (`_laplace_block_lowest`); made dense, it was the unit's width cubed, hours
# for a study of thousands. The dense decompositions are the oracle
# (`ctsem_set_block_eigen!(false)`).
# (`.nested_panel_data()`, `.nested_panel_model()`: helper-nested.R)

test_that("a nested unit's mode search reads its blocks, as its dense oracle does", {
  skip_without_julia()
  withr::local_options(ctsem.laplace.hybrid = FALSE)
  m <- .nested_panel_model("intercepts")
  spec <- suppressMessages(ctFit(.nested_panel_data(), m, backend = "julia",
    intoverpop = "laplace", fit = FALSE, cores = 2))
  expect_null(spec$laplace$hybrid)
  obj <- ctsem:::.ctJuliaObjective(spec, layer = TRUE)
  dims <- ctsem:::.ctJuliaEval("o -> collect(Int, o.units.dims)")
  expect_equal(as.integer(dims(obj)), c(52L, 52L))
  block <- function(on) invisible(ctsem:::.ctJuliaCall(
    "ContinuousTimeSEM.ctsem_set_block_eigen!", on))
  on.exit(block(TRUE), add = TRUE)
  # At two threads, so the unit's leaf blocks are formed at once, which the
  # kept columns once raced on (`_tan_keep_begin!`).
  set.seed(6)
  for (sdx in c(0, 0.3)) {
    x <- rnorm(ctsem:::.ctBackendNpar(spec), 0, sdx)
    block(TRUE)
    a <- ctsem:::.ctJuliaEvaluate(spec, x, layer = TRUE, cores = 2)
    block(FALSE)
    b <- ctsem:::.ctJuliaEvaluate(spec, x, layer = TRUE, cores = 2)
    expect_equal(a$value, b$value, tolerance = 1e-8)
    expect_equal(a$gradient, b$gradient, tolerance = 1e-6)
  }
  # A subset of its units, as the random-effects check takes, for a model with
  # as many units as levels (`affine` and `inert` are per level).
  subset <- ctsem:::.ctJuliaEval(paste0("o -> ContinuousTimeSEM._ctsem_subset_objective(",
    "o, [1]) isa ContinuousTimeSEM.CTSEMLaplaceObjective"))
  expect_true(isTRUE(subset(obj)))
})

# The reflect fixture's people nested in two studies, each study with its own
# random intercept: a study is one unit, searched a block at a time, each
# person's starts given the study's effect (`_laplace_probe_nested`), and a
# better maximum one person's search finds is where the next person's starts
# leave from. Each person is checked against their own maxima given the rest of
# the unit's mode, by a grid refined by Newton and scored by the unit's term.
# The case under test: a person whose better maximum lies past the origin's,
# and one whose better maximum becomes the better one only once an earlier
# person's has moved the study's intercept (person 2 of study 1, by 0.07).
# Persons with sign-ambiguous loadings in two studies of six, each study one
# nested unit: a study intercept over the persons' loadings.
.modes_nested_spec <- function() {
  set.seed(3)
  n <- 12; nt <- 10
  d <- data.frame(id = rep(1:n, each = nt), time = rep(0:(nt - 1), n))
  d$study <- ifelse(d$id <= n / 2, 1, 2)
  eta <- unlist(lapply(1:n, function(i) {
    e <- numeric(nt); e[1] <- 1 + rnorm(1)
    for (t in 2:nt) e[t] <- 0.6 * e[t - 1] + rnorm(1, 0, 0.8)
    e
  }))
  d$Y1 <- -1.2 * eta + rnorm(n * nt, 0, 0.3)
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", id = c("id", "study"),
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix("lambda"),
    DRIFT = matrix(-0.5), DIFFUSION = matrix(0.8), MANIFESTVAR = matrix(0.3),
    MANIFESTMEANS = matrix("mm"), T0MEANS = matrix(1), T0VAR = matrix(1),
    CINT = matrix(0), tipredDefault = FALSE)))
  m$pars$indvarying <- ifelse(m$pars$param %in% "lambda", "id",
    ifelse(m$pars$param %in% "mm", "study", FALSE))
  m$pars$transform[m$pars$param %in% "lambda"] <- "param"
  suppressMessages(ctFit(d, m, backend = "julia", intoverpop = "laplace",
    fit = FALSE, cores = 1))
}

test_that("a nested unit's persons each reach their best maximum given the rest", {
  skip_without_julia()
  n <- 12
  spec <- .modes_nested_spec()
  obj <- ctsem:::.ctJuliaObjective(spec)
  dims <- ctsem:::.ctJuliaEval("o -> collect(Int, o.units.dims)")
  expect_equal(as.integer(dims(obj)), c(7L, 7L))
  x <- numeric(ctsem:::.ctBackendNpar(spec))
  # Scored alone, which maximum the search reached.
  invisible(ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_mode_sum!", FALSE))
  on.exit(invisible(ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_mode_sum!",
    TRUE)), add = TRUE)
  ctsem:::.ctJuliaEvaluate(spec, x, gradient = FALSE, cores = 1L)
  ctsem:::.ctJuliaEval(paste(
    "function _ctsem_test_nested_best(o, theta)",
    "  C = ContinuousTimeSEM",
    "  Ls = C._laplace_popchols(theta, o.spec)",
    "  aws = C._laplace_workspace!(o, Float64, length(theta))",
    "  grid = collect(-6.0:0.01:6.0)",
    "  out = zeros(0, 4)",
    "  for U in eachindex(o.units.members)",
    "    u0 = copy(o.modes[U])",
    "    origin = C._laplace_newton_unit_mode(o, U, theta, Ls, aws, zeros(length(u0))).u",
    "    for blk in o.units.blocks[U]",
    "      isempty(blk.ancestors) && continue",
    "      j = blk.offset + 1",
    "      g = map(grid) do v",
    "        u = copy(u0); u[j] = v",
    "        C._laplace_unit_objective_gradient(o, U, theta, Ls, u, aws).value",
    "      end",
    "      best = -Inf; at = NaN",
    "      for i in 2:length(grid)-1",
    "        (g[i] > g[i-1] && g[i] >= g[i+1]) || continue",
    "        u = copy(u0); u[j] = grid[i]",
    "        r = C._laplace_newton_unit_mode(o, U, theta, Ls, aws, u)",
    "        t = C._laplace_unit_primal(o, U, theta, Ls, copy(r.u), aws).term",
    "        t > best && (best = t; at = r.u[j])",
    "      end",
    "      out = vcat(out, [U u0[j] at origin[j]])",
    "    end",
    "  end",
    "  out",
    "end", sep = "\n"))
  scan <- matrix(as.numeric(ctsem:::.ctJuliaCall("Main._ctsem_test_nested_best", obj,
    ctsem:::.ctJuliaNumericVector(x))), ncol = 4,
    dimnames = list(NULL, c("unit", "found", "best", "origin")))
  expect_equal(nrow(scan), n)
  expect_equal(scan[, "found"], scan[, "best"], tolerance = 1e-3)
  # Exercised: some person's best maximum is not on the side Newton from the
  # origin reaches.
  expect_true(any(sign(scan[, "found"]) != sign(scan[, "origin"])))
})

# Persons with their own random drift under a random study drift, on the layer:
# with wide population sds, short records leave a unit's curvature soft, and the
# soft-direction quadrature's node floor binds. Where it does, the rule and its
# reverse take it from the unit's blocks (the floored eigenpairs' low-rank
# part, `_soft_floored_node_reverse!`), and the dense rule and reverse, which
# take every eigenpair, are their oracle. The point is fixed; the test first
# checks that a node there holds two floored pairs, so it fails if it stops
# exercising them.
test_that("a nested unit's persons combine their own maxima, continuously", {
  skip_without_julia()
  spec <- .modes_nested_spec()
  obj <- ctsem:::.ctJuliaObjective(spec)
  np <- ctsem:::.ctBackendNpar(spec)
  ctsem:::.ctJuliaEval(paste(
    "_ctsem_test_bs_start(o, x) = (ContinuousTimeSEM._ctsem_anchor_start!(o, collect(Float64, x)); nothing)",
    "function _ctsem_test_bs_step(o, x)",
    "  r = ContinuousTimeSEM.ctsem_laplace_evaluate(o, collect(Float64, x); gradient=true)",
    "  ContinuousTimeSEM._ctsem_anchor_here!(o)",
    "  vcat(r.value, r.gradient, sum(length, o.block_modes), sum(length, o.alt_modes))",
    "end",
    "_ctsem_test_bs_value(o, x) = ContinuousTimeSEM.ctsem_laplace_evaluate(o, collect(Float64, x); gradient=false).value",
    "_ctsem_test_bs_clear(o) = (ContinuousTimeSEM.ctsem_laplace_clear_anchors!(o); nothing)",
    sep = "\n"))
  on.exit({
    invisible(ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_mode_sum!", TRUE))
    invisible(ctsem:::.ctJuliaCall("Main._ctsem_test_bs_clear", obj))
  }, add = TRUE)
  jv <- function(f, x) as.numeric(ctsem:::.ctJuliaCall(f, obj,
    ctsem:::.ctJuliaNumericVector(x)))
  # Where every person's own maximum is the best at full weight, the term is
  # the unit's at its mode.
  x <- c(0.3, 0.5, 0, 0)[seq_len(np)]
  on <- ctsem:::.ctJuliaEvaluate(spec, x, gradient = FALSE, cores = 1L)$value
  invisible(ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_mode_sum!", FALSE))
  off <- ctsem:::.ctJuliaEvaluate(spec, x, gradient = FALSE, cores = 1L)$value
  invisible(ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_mode_sum!", TRUE))
  expect_equal(on, off, tolerance = 1e-10)
  # At the population mean every person's other maximum is the person's (no
  # maximum of the whole study is kept), and one of them carries the term.
  x <- numeric(np)
  invisible(jv("Main._ctsem_test_bs_start", x))
  r <- jv("Main._ctsem_test_bs_step", x)
  expect_gt(r[np + 2], 0)
  expect_equal(r[np + 3], 0)
  invisible(ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_mode_sum!", FALSE))
  single <- jv("Main._ctsem_test_bs_value", x)
  invisible(ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_mode_sum!", TRUE))
  expect_gt(r[1] - single, 1e-3)
  # Held, as a run evaluates, the gradient is the value's.
  h <- 1e-6
  fd <- sapply(seq_len(np), function(j) {
    e <- numeric(np); e[j] <- h
    (jv("Main._ctsem_test_bs_value", x + e) - jv("Main._ctsem_test_bs_value", x - e)) /
      (2 * h)
  })
  expect_equal(r[2:(np + 1)], fd, tolerance = 1e-5)
  invisible(ctsem:::.ctJuliaCall("Main._ctsem_test_bs_clear", obj))
  # Walked through the stretch where persons' preferences turn, the anchors
  # moved at every step as an optimiser's accepted steps move them, each
  # step's change against the trapezoid of its end gradients: creases where a
  # person's alternative overtakes its own maximum, which shrink with the
  # step, and no jump, which would not. A study-wide maximum kept beside the
  # persons' credited their flips twice and dropped the term 3 nats here.
  walk <- function(K) {
    x0 <- c(0.15, 0.09375, 0, 0)[seq_len(np)]
    dir <- (c(0.10, 0.1125, 0, 0)[seq_len(np)] - x0) / K
    invisible(jv("Main._ctsem_test_bs_start", x0))
    out <- t(sapply(0:K, function(k) jv("Main._ctsem_test_bs_step", x0 + k * dir)))
    invisible(ctsem:::.ctJuliaCall("Main._ctsem_test_bs_clear", obj))
    g <- out[, 2:(np + 1)] %*% dir
    max(abs(diff(out[, 1]) - (head(g, -1) + tail(g, -1)) / 2))
  }
  coarse <- walk(50)
  fine <- walk(200)
  expect_lt(coarse, 0.05)
  expect_lt(fine, coarse / 2)
})

test_that("a nested unit's soft quadrature takes a binding floor from its blocks", {
  skip_without_julia()
  withr::local_options(ctsem.laplace.hybrid = FALSE)
  m <- .nested_panel_model("drift")
  spec <- suppressMessages(ctFit(.nested_panel_data(3, 8, 12), m, backend = "julia",
    intoverpop = "laplace", fit = FALSE, cores = 2))
  obj <- ctsem:::.ctJuliaObjective(spec, layer = TRUE)
  x <- rep(0, ctsem:::.ctBackendNpar(spec))
  for (lv in spec$laplace$levels) x[lv$sd_index[lv$sd_index > 0]] <- 2.5
  x <- x + 1.2 * as.numeric(ctsem:::.ctJuliaEval(
    "n -> randn(ContinuousTimeSEM.Random.Xoshiro(37), n)")(length(x)))
  # The modes Newton from the origin reaches, the probes off: they find
  # better, stiffer maxima of these units, and the rule is tested where the
  # floor binds, at whichever mode that is.
  compare <- ctsem:::.ctJuliaEval(paste(
    "(o, x) -> begin",
    "  C = ContinuousTimeSEM",
    "  was = C.ctsem_set_mode_probes!(false)",
    "  try C.ctsem_laplace_evaluate(o, x; gradient = false) finally C.ctsem_set_mode_probes!(was) end",
    "  theta = C._laplace_extend(o, collect(Float64, x))",
    "  Ls = C._laplace_popchols(theta, o.spec)",
    "  aws = C._laplace_workspace!(o, Float64, length(theta))",
    "  top = C._laplace_soft_node_curvature_floor_top()",
    "  out = Float64[]",
    "  for U in eachindex(o.units.members)",
    "    u = o.modes[U]; blocks = o.units.blocks[U]",
    "    curvf = z -> C._laplace_unit_curvature(o, U, theta, Ls, z)",
    "    M = C._laplace_copy_blockmatrix(curvf(u))",
    "    F0 = C._soft_shifted_factor(M, blocks, 0.0)",
    "    F0 === nothing && continue",
    "    eig = C._soft_smallest_eigenpair(M, blocks, F0)",
    "    (eig === nothing || eig.lambda >= o.soft_blend_hi) && continue",
    "    most = 0",
    "    for xg in C._gauss_hermite(3)[1]",
    "      abs(xg) < 1e-12 && continue",
    "      ub = u .+ (sqrt(2) * xg / sqrt(max(eig.lambda, 1.0))) .* eig.v",
    "      low = C._soft_node_low_pairs(curvf(ub), blocks, eig.v, top)",
    "      most = max(most, low === nothing ? -1 : length(low.values))",
    "    end",
    "    most > 0 || continue",
    "    inner = C._laplace_unit_objective_gradient(o, U, theta, Ls, u, aws)",
    "    ok, logdetM, _, _ = C._laplace_block_factor(M, blocks)",
    "    fac = C._laplace_soft_quadrature_factored(o, M, blocks, logdetM + 2 * inner.hlogdet,",
    "      inner.value, u, z -> C._laplace_unit_objective_gradient(o, U, theta, Ls, z, aws), curvf)",
    "    prev = C.ctsem_set_soft_quadrature_factored!(false)",
    "    dense = try C._laplace_soft_quadrature_term(o, U, theta, Ls, u, aws) finally",
    "      C.ctsem_set_soft_quadrature_factored!(prev) end",
    "    st = C._laplace_soft_quadrature_setup(o, U, theta, u)",
    "    gf = C._laplace_soft_quadrature_reverse_factored(o, st.z0, st.objf, st.curvb,",
    "      st.blocks, st.G, length(theta))",
    "    gd = C._laplace_soft_quadrature_reverse(o, st.z0, st.objf, st.curvf, st.G, length(theta))",
    "    push!(out, most, fac === nothing ? NaN : fac.value, dense.value,",
    "      gf === nothing ? NaN : maximum(abs.(gf .- gd)), maximum(abs, gd))",
    "  end",
    "  out",
    "end", sep = "\n"))
  r <- matrix(as.numeric(compare(obj, x)), nrow = 5)
  # most floored pairs at a node, block value, dense value, gradient gap, scale
  expect_gt(ncol(r), 0)
  expect_gte(max(r[1, ]), 2)
  expect_equal(r[2, ], r[3, ], tolerance = 1e-10)
  expect_true(all(r[4, ] <= 1e-7 * pmax(1, r[5, ])))
})

# Known modes (`spec$laplace$known`, `_laplace_known_modes`): starts a unit's
# solve takes besides the origin and the probes, kept where a fit saw the
# unit's term jump (`.ctJuliaOptimise()`). The sign fixture at a population
# mean of zero, where every origin is the stationary point between a person's
# two maxima: with the probes off Newton stops there, and given the modes the
# probes find it reaches them without the probes. Nearby the objective with
# them is smooth, its gradient the value's.
test_that("a unit's solve also starts from its known modes", {
  skip_without_julia()
  spec <- suppressMessages(ctFit(.modes_sign_data(), .modes_sign_model(),
    backend = "julia", intoverpop = "laplace", fit = FALSE))
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  probes <- function(on) invisible(ctsem:::.ctJuliaCall(
    "ContinuousTimeSEM.ctsem_set_mode_probes!", on))
  on.exit(probes(TRUE), add = TRUE)
  # Which maximum a known start reaches, scored alone: given one of a
  # person's two, the combination (`_laplace_mode_sum`) would have only that one.
  invisible(ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_mode_sum!", FALSE))
  on.exit(invisible(ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_mode_sum!",
    TRUE)), add = TRUE)
  evaluate <- function(h, x) ctJuliaEvaluate(h, x, gradient = TRUE,
    contributions = TRUE, cores = 1L)
  modes <- ctsem:::.ctJuliaEval("o -> reduce(hcat, o.modes)")
  x <- c(0, 0.5)
  probed <- evaluate(handle, x)
  found <- as.matrix(modes(ctsem:::.ctJuliaObjective(handle)))
  n <- ncol(found)
  known <- spec
  known$laplace$known <- list(units = seq_len(n), values = found)
  kept <- structure(known, class = c("ctJuliaModel", "ctFitModel"))
  # The engine holds what it was sent.
  held <- ctsem:::.ctJuliaEval("o -> sum(length, o.known)")
  expect_equal(held(ctsem:::.ctJuliaObjective(kept)), n)
  probes(FALSE)
  plain <- evaluate(handle, x)
  withknown <- evaluate(kept, x)
  # Without the probes every solve stops at the origin, scored by the soft-direction
  # rule (the first test above), and lands elsewhere than the probes' maxima.
  expect_gt(max(abs(probed$subject_loglik - plain$subject_loglik)), 0.1)
  expect_equal(as.numeric(withknown$subject_loglik),
    as.numeric(probed$subject_loglik), tolerance = 1e-8)
  # Off the stationary point, still without the probes.
  x <- c(0.05, 0.5)
  r <- evaluate(kept, x)
  set.seed(1)
  v <- stats::rnorm(2)
  v <- v / sqrt(sum(v^2))
  h <- 1e-5
  fd <- (evaluate(kept, x + h * v)$value - evaluate(kept, x - h * v)$value) / (2 * h)
  expect_equal(fd, sum(as.numeric(r$gradient) * v), tolerance = 1e-5)
})

# What a jump adds to the known modes (`.ctJuliaKnownModes()`): the columns of
# the units it named, padded with NaN to the longer matrix, none already held.
test_that("a jump's modes join the known ones once each", {
  result <- list(stop_reason = "jump", jump_units = c(3L, 5L),
    jump_point_units = c(3L, 3L, 4L, 5L),
    jump_point_values = matrix(c(1, 2, 3, 4, 5, 6, 7, 8), 2))
  first <- ctsem:::.ctJuliaKnownModes(NULL, result)
  expect_equal(first$units, c(3L, 3L, 5L))
  expect_equal(first$values, matrix(c(1, 2, 3, 4, 7, 8), 2))
  # Seen again, nothing is new.
  expect_null(ctsem:::.ctJuliaKnownModes(first, result))
  # A longer unit's mode pads the rows of those held.
  wider <- list(stop_reason = "jump", jump_units = 7L, jump_point_units = 7L,
    jump_point_values = matrix(c(9, 9, 9), 3))
  both <- ctsem:::.ctJuliaKnownModes(first, wider)
  expect_equal(both$units, c(3L, 3L, 5L, 7L))
  expect_equal(both$values, cbind(rbind(first$values, NaN), c(9, 9, 9)))
  # A mode is its entries before the padding: the same entries padded are the
  # same mode, and a shorter one with the same start is another.
  expect_null(ctsem:::.ctJuliaKnownAdd(both, 3L, matrix(c(1, 2, NaN))))
  expect_equal(ctsem:::.ctJuliaKnownAdd(both, 7L, matrix(c(9, 9, NaN)))$units,
    c(3L, 3L, 5L, 7L, 7L))
  expect_null(ctsem:::.ctJuliaKnownModes(first, list(stop_reason = "jump",
    jump_units = 0L, jump_point_units = 0L, jump_point_values = matrix(0))))
})

# A known mode is used only where its length is the unit's on the objective
# being built: a hybrid's mode, shorter, is passed over by the layer, where
# read as the unit's first entries it once stopped the build.
test_that("a known mode of another length is passed over", {
  skip_without_julia()
  spec <- suppressMessages(ctFit(.modes_sign_data(), .modes_sign_model(),
    backend = "julia", intoverpop = "laplace", fit = FALSE))
  spec$laplace$known <- list(units = c(1L, 1L, 2L),
    values = matrix(c(0.5, NaN, 0.5, 0.7, 0.9, NaN), 2))
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  objective <- ctsem:::.ctJuliaObjective(handle)
  held <- ctsem:::.ctJuliaEval("o -> Int[length(k) for k in o.known]")
  expect_equal(as.integer(held(objective))[1:3], c(1L, 1L, 0L))
})

# Anchors (`_ctsem_anchor_start!`): within an optimisation each unit's solve
# starts from its mode at the last accepted point instead of searching. Near
# that point it reaches the mode a search reaches, so values and gradients
# agree; a fit anchored reaches the optimum a fit searching at every evaluation
# reaches (`ctsem_set_mode_anchors!(false)`, the oracle), and leaves no anchors
# on the objective R caches. The reflect fixture, whose people include one a
# search finds only past the origin.
test_that("anchored solves reach the searched modes, and a fit leaves none", {
  skip_without_julia()
  m <- .modes_sign_model()
  m$pars$value[m$pars$matrix %in% "T0MEANS"] <- 1
  data <- .modes_reflect_data()
  spec <- suppressMessages(ctFit(data, m, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  objective <- ctsem:::.ctJuliaObjective(handle)
  anchor <- ctsem:::.ctJuliaEval(
    "(o, x) -> ContinuousTimeSEM._ctsem_anchor_start!(o, x)")
  clear <- ctsem:::.ctJuliaEval(
    "o -> ContinuousTimeSEM.ctsem_laplace_clear_anchors!(o)")
  held <- ctsem:::.ctJuliaEval("o -> length(o.anchors)")
  x <- c(0.5, 0.5)
  anchor(objective, ctsem:::.ctJuliaNumericVector(x))
  expect_equal(held(objective), length(unique(data$id)))
  near <- x + c(0.01, -0.02)
  a <- ctJuliaEvaluate(handle, near, gradient = TRUE, contributions = TRUE,
    cores = 1L)
  clear(objective)
  expect_equal(held(objective), 0)
  s <- ctJuliaEvaluate(handle, near, gradient = TRUE, contributions = TRUE,
    cores = 1L)
  expect_equal(as.numeric(a$subject_loglik), as.numeric(s$subject_loglik),
    tolerance = 1e-8)
  expect_equal(as.numeric(a$gradient), as.numeric(s$gradient), tolerance = 1e-6)

  fitwith <- function(on) {
    previous <- ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_mode_anchors!", on)
    on.exit(ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_mode_anchors!",
      previous), add = TRUE)
    set.seed(5)
    suppressWarnings(suppressMessages(ctFit(data, m, backend = "julia",
      intoverpop = "laplace", cores = 1L,
      optimcontrol = list(checks = FALSE))))
  }
  anchored <- fitwith(TRUE)
  searched <- fitwith(FALSE)
  expect_equal(anchored$estimate$loglik, searched$estimate$loglik,
    tolerance = 1e-6)
  # The population mean; the sd collapses here, along a flat direction where
  # a fit stops anywhere (raw -4.78 and -4.76 at one setting).
  expect_equal(anchored$estimate$raw[1], searched$estimate$raw[1], tolerance = 1e-4)
  expect_equal(held(ctsem:::.ctJuliaObjective(anchored)), 0)
})

test_that("a checkpoint keeps the point and the anchored modes, and a fit resumes from it", {
  skip_without_julia()
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  data <- .modes_sign_data()
  fitto <- function(model, maxiter) {
    set.seed(5)
    suppressWarnings(ctFit(data, model, backend = "julia", intoverpop = "laplace",
      cores = 1L, optimcontrol = list(checks = FALSE, estonly = TRUE,
        maxiter = maxiter, checkpoint = path)))
  }
  first <- suppressMessages(fitto(.modes_sign_model(), 4L))
  expect_true(file.exists(path))
  saved <- readRDS(path)
  expect_identical(saved$format, "ctsem checkpoint")
  expect_equal(saved$raw, as.numeric(first$estimate$raw))
  # Every unit's anchored mode, each a loading effect, and the mirror maxima
  # the sum over a unit's maxima follows beside them.
  expect_setequal(unique(saved$modes$units), seq_len(length(unique(data$id))))
  expect_true(all(is.finite(saved$modes$values[1, ])))

  # The same call again starts where the first stopped, its search starting
  # from the saved modes too, which it keeps no further.
  expect_message(second <- fitto(.modes_sign_model(), 40L),
    "Resuming from the checkpoint")
  expect_gte(second$estimate$loglik, first$estimate$loglik - 1e-8)
  known <- second$model_spec$laplace$known
  expect_true(is.null(known) || length(known$units) < length(saved$modes$units))

  # Written for another model, the file is not overwritten.
  other <- suppressWarnings(suppressMessages(ctModel(type = "ct", id = "id",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix("lambda"),
    DRIFT = matrix("drift"), DIFFUSION = matrix(0.8), MANIFESTVAR = matrix(0.3),
    MANIFESTMEANS = matrix(0), T0MEANS = matrix(0), T0VAR = matrix(1),
    CINT = matrix(0))))
  other$pars$indvarying <- other$pars$param %in% "lambda"
  expect_error(suppressMessages(fitto(other, 4L)), "written for another model")
  expect_identical(readRDS(path)$npar, saved$npar)
})

# The exact Hessian a finish forms is kept beside the checkpoint, and a resume
# at the same point walks on it rather than forming it again: on the affect
# pilot a run lost after forming its certification's Hessian formed it again
# on resuming, 8 hours at the same point.
test_that("a checkpoint keeps the exact Hessian, and a resume does not form it again", {
  skip_without_julia()
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(c(path, paste0(path, ".hessian"))), add = TRUE)
  data <- .modes_sign_data()
  fitto <- function() {
    set.seed(5)
    suppressWarnings(ctFit(data, .modes_sign_model(), backend = "julia",
      intoverpop = "laplace", cores = 1L,
      optimcontrol = list(checks = FALSE, checkpoint = path)))
  }
  first <- suppressMessages(fitto())
  expect_gte(first$optim$newton_hessians, 1L)
  saved <- readRDS(path)
  kept <- readRDS(paste0(path, ".hessian"))
  expect_identical(kept$id, saved$hessian_id)
  npar <- length(first$estimate$raw)
  expect_equal(dim(kept$hessian), c(npar, npar))
  expect_identical(kept$state, saved$state)

  expect_message(second <- fitto(), "with the exact Hessian it kept")
  expect_equal(second$optim$newton_hessians, 0L)
  expect_equal(second$estimate$loglik, first$estimate$loglik, tolerance = 1e-6)
  expect_true(isTRUE(second$optim$converged))
  # Still named by the checkpoint the resume wrote, for a resume after it.
  expect_identical(readRDS(path)$hessian_id, kept$id)
})

# After the fit, what reads it holds the modes its run ended with
# (R/ctJuliaFitModes.R), and the fit keeps them, found at its estimate.
test_that("a fit keeps the modes its run ended with", {
  skip_without_julia()
  set.seed(5)
  fit <- suppressWarnings(suppressMessages(ctFit(.modes_sign_data(),
    .modes_sign_model(), backend = "julia", intoverpop = "laplace", cores = 1L,
    optimcontrol = list(checks = FALSE))))
  modes <- fit$laplace$modes
  expect_setequal(unique(modes$units), seq_len(length(unique(.modes_sign_data()$id))))
  expect_equal(ncol(modes$values), length(modes$units))
  expect_equal(modes$at, as.numeric(fit$estimate$raw))
})

# Held, the objective is a function of the parameters given the modes: every
# evaluation starts each unit from them and follows each maximum from where
# they were found, whatever was evaluated before, and no search from the origin
# adds a maximum they lack. At a population mean of 0.3 each person on the sign
# fixture has two maxima.
test_that("held modes give one function of the parameters, and holds nest", {
  skip_without_julia()
  spec <- suppressMessages(ctFit(.modes_sign_data(), .modes_sign_model(),
    backend = "julia", intoverpop = "laplace", fit = FALSE))
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  objective <- ctsem:::.ctJuliaObjective(handle)
  module <- ctsem:::.ctJuliaModule(spec$project)
  value <- function(x) as.numeric(ctsem:::.ctJuliaGet(module$ctsem_laplace_evaluate(
    objective, ctsem:::.ctJuliaNumericVector(x), gradient = FALSE))$value)
  ctsem:::.ctJuliaEval(paste(
    "function _ctsem_test_held_modes(o, x)",
    "  C = ContinuousTimeSEM",
    "  C.ctsem_laplace_clear_anchors!(o)",
    "  C.ctsem_laplace_evaluate(o, x; gradient = false)",
    "  units = Int[]; cols = Vector{Float64}[]",
    "  for U in eachindex(o.modes)",
    "    push!(units, U); push!(cols, copy(o.modes[U]))",
    "    for m in o.alt_modes[U]; push!(units, U); push!(cols, copy(m)); end",
    "  end",
    "  V = fill(NaN, maximum(length, cols), length(cols))",
    "  for j in eachindex(cols); V[1:length(cols[j]), j] .= cols[j]; end",
    "  (units = units, values = V)",
    "end", sep = "\n"))
  x0 <- c(0.3, 0.5)
  r <- ctsem:::.ctJuliaGet(ctsem:::.ctJuliaCall("Main._ctsem_test_held_modes",
    objective, ctsem:::.ctJuliaNumericVector(x0)))
  modes <- list(units = as.integer(r$units), values = as.matrix(r$values), at = x0)
  expect_gt(length(modes$units), length(unique(modes$units)))
  points <- list(x0, x0 + c(0.2, 0.1), x0 - c(0.2, 0.1))
  # Each point's value under a hold, released when this returns.
  held <- function(modes) {
    expect_equal(ctsem:::.ctJuliaHoldModes(spec, modes), 12L)
    vapply(points, value, numeric(1))
  }
  searched <- vapply(points, value, numeric(1))
  full <- held(modes)
  expect_equal(full[1:2], searched[1:2], tolerance = 1e-8)
  # One function: a point's value whatever was evaluated before it.
  again <- local({
    ctsem:::.ctJuliaHoldModes(spec, modes)
    c(value(points[[3]]), value(points[[2]]))
  })
  expect_identical(again[2], full[2])
  # Released, the search again.
  expect_identical(value(points[[2]]), searched[2])
  # Each unit held at its lower maximum alone: the hold is what it holds, and
  # nothing searches for the higher.
  lower <- !duplicated(modes$units, fromLast = TRUE)
  some <- list(units = modes$units[lower],
    values = modes$values[, lower, drop = FALSE], at = x0)
  short <- held(some)
  expect_true(all(short < searched - 1))
  # A hold inside another is released alone, and the one under it holds again
  # with its own modes.
  nested <- local({
    ctsem:::.ctJuliaHoldModes(spec, modes)
    inner <- local({
      ctsem:::.ctJuliaHoldModes(spec, some)
      value(points[[1]])
    })
    c(inner, value(points[[1]]))
  })
  expect_equal(nested, c(short[1], full[1]), tolerance = 1e-8)
})
