# Own variances on a reduced-rank level: `ctFit(..., popunique =)`.
#
# A level's covariance under `poprank = r` with own variances is
# `L L' + diag(psi^2)`, the factor-analytic form; `poprank = 0` with every
# effect selected is a diagonal covariance. The engine builds it as one root,
# `F = [L | Psi]` (`_laplace_poploading`), each own variance transformed as a
# full-rank sd. These tests check the covariance against one computed here from
# the raw values, the diagonal case against a full-rank level whose
# correlations are stated zero (the same model by another route), the
# gradient against forward mode, the names, and the refusals.

.pu_model <- function(levels = "id") {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    id = levels, manifestNames = c("Y1", "Y2", "Y3"), latentNames = "eta1",
    LAMBDA = matrix(c(1, 0.8, 1.2), 3, 1), DRIFT = matrix(-0.5),
    DIFFUSION = matrix(0.5), MANIFESTVAR = diag(0.4, 3),
    MANIFESTMEANS = matrix(c("mm1", "mm2", "mm3"), 3, 1),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0))))
  m$pars$indvarying <- m$pars$param %in% c("mm1", "mm2", "mm3")
  if ("study" %in% levels) m$pars$indvarying_study <- m$pars$param %in% c("mm1", "mm2")
  m
}

.pu_data <- function(study = FALSE) {
  set.seed(11)
  n <- 24; t <- 5
  d <- data.frame(id = rep(seq_len(n), each = t), time = rep(seq_len(t) - 1, n))
  person <- matrix(rnorm(n * 3), n, 3)
  for (j in 1:3) d[[paste0("Y", j)]] <- person[d$id, j] + rnorm(n * t)
  if (study) d$study <- ceiling(d$id / 6)
  d
}

.pu_spec <- function(d, m, ...) suppressMessages(ctFit(d, m, backend = "julia",
  intoverpop = "laplace", fit = FALSE, ...))

.pu_names <- function(spec) {
  npar <- ctsem:::.ctBackendNpar(spec)
  ctsem:::.ctBackendRawParameterNames(list(model_spec = spec), npar)
}

.pu_eval <- function(spec, x, method = "adjoint") ctJuliaEvaluate(
  structure(spec, class = c("ctJuliaModel", "ctFitModel")), x,
  gradient = TRUE, gradient_method = method, cores = 1L)

test_that("poprank 0 with own variances is the diagonal covariance", {
  skip_without_julia()
  d <- .pu_data()
  diag_spec <- .pu_spec(d, .pu_model(), poprank = 0, popunique = TRUE)
  # The same model the long way: full rank, every correlation coordinate stated 0.
  m <- ctsem:::.ctModelRawPopVarSync(.pu_model())
  for (a in c("mm1", "mm2", "mm3")) for (b in c("mm1", "mm2", "mm3")) {
    if (a != b) m$RAWPOPVAR[a, b] <- 0
  }
  full_spec <- .pu_spec(d, m)
  nd <- .pu_names(diag_spec)
  nf <- .pu_names(full_spec)
  expect_setequal(grep("^popunique_", nd, value = TRUE),
    paste0("popunique_", c("mm1", "mm2", "mm3")))
  expect_length(grep("^poploading_", nd), 0L)
  expect_equal(length(nd), length(nf))
  set.seed(3)
  xf <- rnorm(length(nf), 0, 0.3)
  xd <- xf[match(sub("^popunique_", "popsd_", nd), nf)]
  expect_false(anyNA(xd))
  rf <- .pu_eval(full_spec, xf)
  rd <- .pu_eval(diag_spec, xd)
  expect_equal(rd$value, rf$value, tolerance = 1e-9)
  expect_equal(as.numeric(rd$gradient),
    as.numeric(rf$gradient)[match(sub("^popunique_", "popsd_", nd), nf)],
    tolerance = 1e-7)
})

test_that("rank 1 with own variances builds L L' + diag(psi^2), and its gradient matches forward mode", {
  skip_without_julia()
  d <- .pu_data()
  spec <- .pu_spec(d, .pu_model(), poprank = 1, popunique = "MANIFESTMEANS")
  nm <- .pu_names(spec)
  set.seed(5)
  x <- rnorm(length(nm), 0, 0.3)
  level <- spec$laplace$levels[[1]]
  expect_equal(level$rank, 1L)
  expect_true(all(level$unique_index > 0L))
  # The covariance the engine builds against one written out here from the raw
  # values: loadings times sdscale, own sds through the full-rank sd transform.
  L <- x[level$load_index] * level$sd_scale
  psi <- log1p(exp(2 * x[level$unique_index] - 1)) * level$sd_scale + 1e-10
  expected <- tcrossprod(L) + diag(psi^2)
  module <- ctsem:::.ctJuliaModule(spec$project)
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  got <- ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_popcov(
    ctsem:::.ctJuliaObjective(handle), ctsem:::.ctJuliaNumericVector(x), 1L))
  expect_equal(unname(as.matrix(got)), unname(expected), tolerance = 1e-10)
  adjoint <- .pu_eval(spec, x, "adjoint")
  forward <- .pu_eval(spec, x, "forward")
  expect_equal(adjoint$value, forward$value, tolerance = 1e-10)
  expect_equal(as.numeric(adjoint$gradient), as.numeric(forward$gradient),
    tolerance = 1e-6)
  # Each own variance moves the objective: a coordinate the engine ignored
  # would have a gradient of exactly zero.
  own <- grep("^popunique_", nm)
  expect_true(all(abs(as.numeric(adjoint$gradient)[own]) > 0))
})

test_that("own variances are named per level and selected by matrix or parameter", {
  skip_without_julia()
  d <- .pu_data(study = TRUE)
  spec <- .pu_spec(d, .pu_model(c("id", "study")), poprank = c(id = 1, study = 1),
    popunique = list(id = "MANIFESTMEANS", study = "mm2"))
  nm <- .pu_names(spec)
  expect_setequal(grep("^popunique_", nm, value = TRUE),
    c(paste0("popunique_", c("mm1", "mm2", "mm3"), ".id"), "popunique_mm2.study"))
  # The layout is contiguous: every raw position is named exactly once.
  expect_false(anyNA(nm))
  expect_false(any(duplicated(nm)))
})

test_that("popunique is refused where it cannot apply, by name", {
  skip_without_julia()
  d <- .pu_data()
  m <- .pu_model()
  expect_error(.pu_spec(d, m, popunique = TRUE), "reduced rank")
  expect_error(.pu_spec(d, m, poprank = NA, popunique = TRUE), "reduced rank")
  expect_error(.pu_spec(d, m, poprank = 3, popunique = TRUE), "reduced rank")
  expect_error(.pu_spec(d, m, poprank = 0, popunique = "mm1"), "every parameter")
  expect_error(.pu_spec(d, m, poprank = 1, popunique = "LAMBDA"), "neither a matrix")
  expect_error(.pu_spec(d, m, poprank = 1, popunique = c(nolevel = TRUE)), "no level")
  expect_error(suppressMessages(ctFit(d, m, backend = "julia",
    intoverpop = "augmented", fit = FALSE, poprank = 1, popunique = TRUE)),
    "popunique requires")
  expect_error(suppressMessages(ctFit(d, m, backend = "stan", fit = FALSE,
    popunique = TRUE)), "popunique requires")
  # FALSE asks for nothing.
  expect_silent(nofit <- .pu_spec(d, m, popunique = FALSE))
})

test_that("own variances of items a subject never answers are skipped, exactly", {
  skip_without_julia()
  # Half the subjects never answer Y3, so its own variance cannot reach their
  # data; the engine skips that coordinate for them. Pruning on and off must
  # give the same value and gradient, on a single-level unit and on the block
  # route of a nested one (forced at any unit size by the block threshold).
  d <- .pu_data(study = TRUE)
  d$Y3[d$id > 12] <- NA
  set.seed(7)
  threshold <- NULL
  on.exit({
    ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_set_laplace_prune!(true)")
    if (!is.null(threshold)) ctsem:::.ctJuliaEval(sprintf(
      "ContinuousTimeSEM.ctsem_set_block_threshold!(%d)", threshold))
  }, add = TRUE)
  for (levels in list("id", c("id", "study"))) {
    rank <- if (length(levels) == 1L) 1 else c(id = 1, study = 1)
    pu <- if (length(levels) == 1L) "MANIFESTMEANS" else list(id = "MANIFESTMEANS")
    spec <- .pu_spec(d, .pu_model(levels), poprank = rank, popunique = pu)
    inert <- spec$laplace$inert
    expect_false(is.null(inert))
    # Y3's intercept (effect 3), for subjects 13-24 at the subject level.
    expect_setequal(inert[inert[, "level"] == 1, "group"], 13:24)
    expect_true(all(inert[, "effect"] == 3L))
    x <- rnorm(length(.pu_names(spec)), 0, 0.3)
    prune <- function(on) ctsem:::.ctJuliaEval(sprintf(
      "ContinuousTimeSEM.ctsem_set_laplace_prune!(%s)", tolower(on)))
    if (length(levels) > 1L) {
      threshold <- as.integer(ctsem:::.ctJuliaEval(
        "ContinuousTimeSEM._LAPLACE_BLOCK_THRESHOLD[]"))
      ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_set_block_threshold!(1)")
    }
    prune(FALSE)
    full <- .pu_eval(spec, x)
    prune(TRUE)
    pruned <- .pu_eval(spec, x)
    expect_equal(pruned$value, full$value, tolerance = 1e-10)
    expect_equal(as.numeric(pruned$gradient), as.numeric(full$gradient),
      tolerance = 1e-8)
  }
})

test_that("own intercept variances the filter integrates match the Laplace coordinates", {
  skip_without_julia()
  # A linear Gaussian model with own variances on intercepts: the Laplace layer
  # is exact over them and so is the filter's closed form (`CTSEMHybrid`), so
  # taking them out of the layer, the default, and leaving them in must give
  # the same value and gradient. Missing items exercise the per-row observed
  # subset, rank 0 a level left with no Laplace coordinate at all, and the
  # study level that a nested level is untouched.
  d <- .pu_data(study = TRUE)
  d$Y3[d$id > 12 & d$time > 1] <- NA
  d$Y1[d$id %% 3 == 0 & d$time == 2] <- NA
  old <- options(ctsem.laplace.hybrid = TRUE)
  on.exit(options(old), add = TRUE)
  set.seed(9)
  cases <- list(
    list(levels = "id", rank = 1, pu = "MANIFESTMEANS"),
    list(levels = "id", rank = 0, pu = TRUE),
    list(levels = c("id", "study"), rank = c(id = 1, study = 1),
      pu = list(id = "MANIFESTMEANS")))
  for (case in cases) {
    m <- .pu_model(case$levels)
    options(ctsem.laplace.hybrid = TRUE)
    hy <- .pu_spec(d, m, poprank = case$rank, popunique = case$pu)
    options(ctsem.laplace.hybrid = FALSE)
    la <- .pu_spec(d, m, poprank = case$rank, popunique = case$pu)
    expect_equal(hy$laplace$hybrid$effect, 1:3)
    expect_null(la$laplace$hybrid)
    x <- rnorm(length(.pu_names(hy)), 0, 0.3)
    rh <- .pu_eval(hy, x)
    rl <- .pu_eval(la, x)
    expect_equal(rh$value, rl$value, tolerance = 1e-9)
    expect_equal(as.numeric(rh$gradient), as.numeric(rl$gradient), tolerance = 1e-7)
  }
})

test_that("integrated own intercepts beside a shared effect on the dynamics: gradient matches forward mode", {
  skip_without_julia()
  # The shared effect loads on DRIFT, so the Laplace layer over it is an
  # approximation and the hybrid is not the same objective as the layer over
  # everything; its adjoint must still agree with forward mode through it.
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    id = "id", manifestNames = c("Y1", "Y2", "Y3"), latentNames = "eta1",
    LAMBDA = matrix(c(1, 0.8, 1.2), 3, 1), DRIFT = matrix("drift"),
    DIFFUSION = matrix(0.5), MANIFESTVAR = diag(0.4, 3),
    MANIFESTMEANS = matrix(c("mm1", "mm2", "mm3"), 3, 1),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0))))
  m$pars$indvarying <- m$pars$param %in% c("drift", "mm1", "mm2", "mm3")
  d <- .pu_data()
  d$Y2[d$id > 16] <- NA
  spec <- .pu_spec(d, m, poprank = 1, popunique = "MANIFESTMEANS")
  expect_equal(length(spec$laplace$hybrid$effect), 3L)
  set.seed(13)
  x <- rnorm(length(.pu_names(spec)), 0, 0.3)
  adjoint <- .pu_eval(spec, x, "adjoint")
  forward <- .pu_eval(spec, x, "forward")
  expect_equal(adjoint$value, forward$value, tolerance = 1e-10)
  expect_equal(as.numeric(adjoint$gradient), as.numeric(forward$gradient),
    tolerance = 1e-6)
})

test_that("a fit with own variances converges and reports them", {
  skip_without_julia()
  skip_on_cran()
  d <- .pu_data()
  fit <- suppressMessages(suppressWarnings(ctFit(d, .pu_model(), backend = "julia",
    intoverpop = "laplace", poprank = 1, popunique = TRUE, priors = TRUE,
    cores = 1)))
  expect_true(isTRUE(fit$optim$converged))
  s <- summary(fit)
  expect_equal(nrow(s$popsd), 3L)
  # Each item's person means vary independently in the data, so the own
  # variances carry the spread: no sd collapses to zero.
  expect_true(all(s$popsd[, "mean"] > 0.3))
})
