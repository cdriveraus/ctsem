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
  # With a study level: mm1 and mm2 at both levels, mm3 between subjects only.
  if ("study" %in% levels) m$pars$indvarying <- ifelse(m$pars$param %in%
    c("mm1", "mm2"), TRUE, ifelse(m$pars$param %in% "mm3", levels[1], FALSE))
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

# Shared effects on DRIFT and DIFFUSION beside own intercepts, and data with one
# item missing for a third of the subjects: the deviations' curvature depends
# on the shared effect, so the filter's profile and the layer must agree through
# that dependence.
.pu_dyn_model <- function() {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    id = "id", manifestNames = c("Y1", "Y2", "Y3"), latentNames = "eta1",
    LAMBDA = matrix(c(1, 0.8, 1.2), 3, 1), DRIFT = matrix("drift"),
    DIFFUSION = matrix("diffusion"), MANIFESTVAR = diag(0.4, 3),
    MANIFESTMEANS = matrix(c("mm1", "mm2", "mm3"), 3, 1),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0))))
  m$pars$indvarying <- m$pars$param %in% c("drift", "diffusion", "mm1", "mm2", "mm3")
  m
}

.pu_dyn_data <- function() {
  d <- .pu_data()
  d$Y2[d$id > 16] <- NA
  d
}

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

test_that("an asymptote reaches its own binary item, not the manifest its row number names", {
  skip_without_julia()
  # Fit-free. ASYMPTOTES row 2 is the second binary indicator, Y3; read as a
  # manifest number it would be Y2. With Y2 missing for subjects 13-24 and Y3
  # for 1-6, Y3's guessing parameter is inert for 1-6 alone.
  set.seed(11)
  n <- 24; t <- 6
  d <- data.frame(id = rep(seq_len(n), each = t), time = rep(seq_len(t) - 1, n))
  eta <- rnorm(n)[d$id] + rnorm(n * t, 0, 0.5)
  d$Y1 <- rbinom(n * t, 1, plogis(eta))
  d$Y2 <- eta + rnorm(n * t, 0, 0.6)
  d$Y3 <- rbinom(n * t, 1, plogis(0.8 * eta))
  d$Y2[d$id > 12] <- NA
  d$Y3[d$id <= 6] <- NA
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", id = "id",
    manifestNames = c("Y1", "Y2", "Y3"), latentNames = "eta1",
    manifesttype = c(1, 0, 1), LAMBDA = matrix(c(1, 0.8, 0.8), 3, 1),
    DRIFT = matrix(-0.5), DIFFUSION = matrix(0.5), MANIFESTVAR = diag(c(0, 0.4, 0)),
    MANIFESTMEANS = matrix(c("mm1", "mm2", "mm3"), 3, 1),
    ASYMPTOTES = rbind(c(NA, NA), c("g3", NA)),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0))))
  m$pars$indvarying <- m$pars$param %in% c("mm2", "g3")
  spec <- .pu_spec(d, m, poprank = 1, popunique = TRUE)
  inert <- spec$laplace$inert
  effects <- spec$laplace$levels[[1]]$param
  groups <- function(name) sort(unname(inert[effects[inert[, "effect"]] == name, "group"]))
  expect_equal(groups("g3"), 1:6)
  expect_equal(groups("mm2"), 13:24)
})

test_that("any measurement effect reaches only its own items, and is skipped where unanswered", {
  skip_without_julia()
  # Y2 and Y3 share an intercept mmA, and Y3 adds mmB times a time dependent
  # predictor, both held in PARS and read by the intercept cells' expressions;
  # Y3 also has a random loading and a random error variance. Everything of
  # Y3's own cannot reach the subjects who never answer it -- reading a
  # predictor joins nothing -- while mmA reaches Y2 as well, which everyone
  # answers. With Y3's error correlated with Y2's, its error variance reaches
  # Y2 through the covariance's construction, and is no longer inert; the
  # loading still is, since the update reads only the observed rows.
  d <- .pu_data()
  d$tp <- rnorm(nrow(d))
  d$Y3[d$id > 12] <- NA
  on.exit(ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_set_laplace_prune!(true)"),
    add = TRUE)
  for (coupled in c(FALSE, TRUE)) {
    m <- suppressWarnings(suppressMessages(ctModel(type = "ct", id = "id",
      manifestNames = c("Y1", "Y2", "Y3"), latentNames = "eta1",
      TDpredNames = "tp", TDPREDEFFECT = matrix(0),
      LAMBDA = matrix(c(1, 0.8, "lam3||TRUE"), 3, 1), DRIFT = matrix(-0.5),
      DIFFUSION = matrix(0.5),
      MANIFESTVAR = matrix(c(0.4, 0, 0, 0, 0.4, if (coupled) "mvc" else 0,
        0, 0, "mv3||TRUE"), 3),
      MANIFESTMEANS = matrix(c("mm1", "mmA", "mmA + mmB * tp"), 3, 1),
      PARS = matrix(c("mmA||TRUE", "mmB||TRUE")),
      T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0))))
    spec <- .pu_spec(d, m, poprank = 1, popunique = TRUE)
    inert <- spec$laplace$inert
    expect_false(is.null(inert))
    expect_setequal(inert[, "group"], 13:24)
    effects <- spec$laplace$levels[[1]]$param
    expect_setequal(effects[unique(inert[, "effect"])],
      c("mmB", "lam3", if (!coupled) "mv3"))
    set.seed(8)
    x <- rnorm(length(.pu_names(spec)), 0, 0.3)
    pruned <- .pu_eval(spec, x)
    ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_set_laplace_prune!(false)")
    full <- .pu_eval(spec, x)
    ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_set_laplace_prune!(true)")
    expect_equal(pruned$value, full$value, tolerance = 1e-10, info = coupled)
    expect_equal(as.numeric(pruned$gradient), as.numeric(full$gradient),
      tolerance = 1e-8, info = coupled)
  }
})

test_that("a process measured by one item reaches that item alone, and is skipped where unanswered", {
  skip_without_julia()
  # eta2 is measured by Y3 alone and fed by eta1, not the reverse, so its
  # drift, intercept, noise, initial mean and predictor effect cannot reach
  # the subjects who never answer Y3 -- linear, and with eta2's noise
  # depending on eta1. eta1's predictor effect reaches everyone. Once eta2
  # feeds what Y1 or Y2 read, nothing of it is inert.
  d <- .pu_data()
  d$tp <- rnorm(nrow(d))
  d$Y3[d$id > 12] <- NA
  build <- function(drift = matrix(c("d11", "d21", 0, "d22"), 2),
    diffusion = matrix(c("q1", 0, 0, "q2"), 2), manifestvar = diag(0.4, 3)) {
    m <- suppressWarnings(suppressMessages(ctModel(type = "ct", id = "id",
      manifestNames = c("Y1", "Y2", "Y3"), latentNames = c("eta1", "eta2"),
      TDpredNames = "tp", TDPREDEFFECT = matrix(c("tpe1", "tpe2"), 2, 1),
      LAMBDA = matrix(c(1, 0.8, 0, 0, 0, 1), 3, 2), DRIFT = drift,
      DIFFUSION = diffusion, MANIFESTVAR = manifestvar,
      MANIFESTMEANS = matrix(c("mm1", "mm2", "mm3"), 3, 1),
      T0VAR = diag(1, 2), T0MEANS = matrix(c("t0m1", "t0m2"), 2),
      CINT = matrix(c(0, "cint2"), 2))))
    m$pars$indvarying <- m$pars$param %in% c("mm1", "d21", "d22", "cint2", "q2",
      "t0m2", "tpe1", "tpe2")
    .pu_spec(d, m, poprank = 1, popunique = TRUE)
  }
  inert_effects <- function(spec) {
    inert <- spec$laplace$inert
    if (is.null(inert)) return(character())
    expect_setequal(inert[, "group"], 13:24)
    spec$laplace$levels[[1]]$param[unique(inert[, "effect"])]
  }
  expect_length(inert_effects(build(drift = matrix(c("d11", "d21", "d12", "d22"), 2))), 0)
  expect_length(inert_effects(build(
    diffusion = matrix(c("0.5 * exp(0.1 * state[2])", 0, 0, "q2"), 2))), 0)
  expect_length(inert_effects(build(manifestvar = matrix(c("0.4 * exp(0.1 * state[2])",
    0, 0, 0, 0.4, 0, 0, 0, 0.4), 3))), 0)

  on.exit(ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_set_laplace_prune!(true)"),
    add = TRUE)
  for (noise in c("q2", "0.5 * exp(0.1 * state[1])")) {
    spec <- build(diffusion = matrix(c("q1", 0, 0, noise), 2))
    expect_setequal(inert_effects(spec),
      c("d21", "d22", "cint2", "t0m2", "tpe2", if (noise == "q2") "q2"))
    set.seed(9)
    x <- rnorm(length(.pu_names(spec)), 0, 0.3)
    pruned <- .pu_eval(spec, x)
    ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_set_laplace_prune!(false)")
    full <- .pu_eval(spec, x)
    ctsem:::.ctJuliaEval("ContinuousTimeSEM.ctsem_set_laplace_prune!(true)")
    expect_equal(pruned$value, full$value, tolerance = 1e-10, info = noise)
    expect_equal(as.numeric(pruned$gradient), as.numeric(full$gradient),
      tolerance = 1e-8, info = noise)
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

test_that("integrated own intercepts beside shared effects on the dynamics are still Laplace exactly", {
  skip_without_julia()
  # The shared effect loads on DRIFT and DIFFUSION, so the deviations'
  # curvature depends on it and the Laplace layer over everything is an
  # approximation. The filter's profile of the deviations, with their logdet
  # added at the mode, is that same approximation: the joint mode, the Schur
  # complement for the curvature, and the joint logdet. So the value and
  # gradient must match the layer holding them as coordinates, and the adjoint
  # must match forward mode through the logdet's dependence on the mode.
  m <- .pu_dyn_model()
  d <- .pu_dyn_data()
  old <- options(ctsem.laplace.hybrid = TRUE)
  on.exit(options(old), add = TRUE)
  spec <- .pu_spec(d, m, poprank = 1, popunique = "MANIFESTMEANS")
  options(ctsem.laplace.hybrid = FALSE)
  layer <- .pu_spec(d, m, poprank = 1, popunique = "MANIFESTMEANS")
  expect_equal(length(spec$laplace$hybrid$effect), 3L)
  set.seed(13)
  for (draw in 1:2) {
    x <- rnorm(length(.pu_names(spec)), 0, 0.3)
    adjoint <- .pu_eval(spec, x, "adjoint")
    forward <- .pu_eval(spec, x, "forward")
    full <- .pu_eval(layer, x, "adjoint")
    expect_equal(adjoint$value, full$value, tolerance = 1e-9)
    expect_equal(as.numeric(adjoint$gradient), as.numeric(full$gradient),
      tolerance = 1e-7)
    expect_equal(adjoint$value, forward$value, tolerance = 1e-10)
    expect_equal(as.numeric(adjoint$gradient), as.numeric(forward$gradient),
      tolerance = 1e-6)
  }
})

test_that("a measurement error written in a time dependent predictor keeps the deviations integrated, exactly", {
  skip_without_julia()
  # A study's known response rounding as a per-row variance q added to each
  # item's error, MANIFESTVAR = sqrt(q + sd^2): a calculation, but in a
  # predictor and parameters only, so each row's covariances and gains are fixed
  # and the filter still integrates the own intercept deviations exactly. The
  # hybrid applies, and equals the layer holding them as coordinates.
  error <- matrix("0", 3, 3)
  diag(error) <- c("sqrt(q + mvar1^2)", "sqrt(q + mvar2^2)", "sqrt(q + mvar3^2)")
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    id = "id", manifestNames = c("Y1", "Y2", "Y3"), latentNames = "eta1",
    TDpredNames = "q", TDPREDEFFECT = matrix(0, 1, 1),
    LAMBDA = matrix(c(1, 0.8, 1.2), 3, 1), DRIFT = matrix("drift"),
    DIFFUSION = matrix("diffusion"),
    PARS = c("mvar1|1e-10 + 5 * log1p_exp(2 * param)",
      "mvar2|1e-10 + 5 * log1p_exp(2 * param)",
      "mvar3|1e-10 + 5 * log1p_exp(2 * param)"),
    MANIFESTVAR = error,
    MANIFESTMEANS = matrix(c("mm1", "mm2", "mm3"), 3, 1),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0))))
  m$pars$indvarying <- m$pars$param %in% c("drift", "diffusion", "mm1", "mm2", "mm3")
  d <- .pu_dyn_data()
  d$q <- ifelse(d$id <= 12, (10 / 6)^2 / 12, 0.1^2 / 12)
  old <- options(ctsem.laplace.hybrid = TRUE)
  on.exit(options(old), add = TRUE)
  spec <- .pu_spec(d, m, poprank = 1, popunique = "MANIFESTMEANS")
  options(ctsem.laplace.hybrid = FALSE)
  layer <- .pu_spec(d, m, poprank = 1, popunique = "MANIFESTMEANS")
  expect_equal(length(spec$laplace$hybrid$effect), 3L)
  set.seed(17)
  for (draw in 1:2) {
    x <- rnorm(length(.pu_names(spec)), 0, 0.3)
    adjoint <- .pu_eval(spec, x, "adjoint")
    forward <- .pu_eval(spec, x, "forward")
    full <- .pu_eval(layer, x, "adjoint")
    expect_equal(adjoint$value, full$value, tolerance = 1e-9)
    expect_equal(as.numeric(adjoint$gradient), as.numeric(full$gradient),
      tolerance = 1e-7)
    expect_equal(as.numeric(adjoint$gradient), as.numeric(forward$gradient),
      tolerance = 1e-6)
  }
  # A calculation reading a state still keeps the deviations in the layer.
  m2 <- m
  m2$pars$param[m2$pars$matrix == "MANIFESTVAR" & m2$pars$row == 1 &
    m2$pars$col == 1] <- "sqrt(q + mvar1^2 + 0.01 * state[1]^2)"
  options(ctsem.laplace.hybrid = TRUE)
  expect_null(.pu_spec(d, m2, poprank = 1, popunique = "MANIFESTMEANS")$laplace$hybrid)

  # The same rounding read from a time independent predictor: one value per
  # subject, so exactly the time dependent version's value and gradient, with
  # the deviations still integrated. The stan backend has no such
  # reference.
  mt <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    id = "id", manifestNames = c("Y1", "Y2", "Y3"), latentNames = "eta1",
    TIpredNames = "r", tipredDefault = FALSE,
    LAMBDA = matrix(c(1, 0.8, 1.2), 3, 1), DRIFT = matrix("drift"),
    DIFFUSION = matrix("diffusion"),
    PARS = c("mvar1|1e-10 + 5 * log1p_exp(2 * param)",
      "mvar2|1e-10 + 5 * log1p_exp(2 * param)",
      "mvar3|1e-10 + 5 * log1p_exp(2 * param)"),
    MANIFESTVAR = {
      e <- matrix("0", 3, 3)
      diag(e) <- c("sqrt(r + mvar1^2)", "sqrt(r + mvar2^2)", "sqrt(r + mvar3^2)")
      e
    },
    MANIFESTMEANS = matrix(c("mm1", "mm2", "mm3"), 3, 1),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0))))
  mt$pars$indvarying <- mt$pars$param %in% c("drift", "diffusion", "mm1", "mm2", "mm3")
  dt <- d
  dt$r <- dt$q
  tispec <- .pu_spec(dt, mt, poprank = 1, popunique = "MANIFESTMEANS")
  expect_equal(length(tispec$laplace$hybrid$effect), 3L)
  expect_equal(.pu_names(tispec), .pu_names(spec))
  set.seed(17)
  for (draw in 1:2) {
    x <- rnorm(length(.pu_names(spec)), 0, 0.3)
    viatd <- .pu_eval(spec, x)
    viati <- .pu_eval(tispec, x)
    expect_equal(viati$value, viatd$value, tolerance = 1e-12)
    expect_equal(as.numeric(viati$gradient), as.numeric(viatd$gradient),
      tolerance = 1e-10)
  }
  # Population matrices read the predictor at zero, or at a value supplied,
  # though it has no effect on any parameter.
  classed <- structure(tispec, class = c("ctJuliaModel", "ctFitModel"))
  at0 <- ctBackendParMatrices(classed, raw = x, trim = FALSE)$MANIFESTVAR
  at5 <- ctBackendParMatrices(classed, raw = x, tipreds = 0.5, trim = FALSE)$MANIFESTVAR
  expect_equal(unname(diag(at5)^2 - diag(at0)^2), rep(0.5, 3), tolerance = 1e-10)
  # A subject's missing value of a predictor a cell reads is imputed, and the
  # cell reads the imputed value; or sampled, a parameter the cell reads, whose
  # gradient comes back through the cell. The model has no effect of the
  # predictor, so without that route only the imputation's own term would.
  gappy <- dt
  gappy$r[gappy$id == 3] <- NA
  imputed <- suppressWarnings(.pu_spec(gappy, mt, poprank = 1,
    popunique = "MANIFESTMEANS"))
  expect_equal(imputed$ti_imputed$subject, 3L)
  expect_true(is.finite(imputed$tipred_data[3, 1]))
  sampled <- suppressWarnings(.pu_spec(gappy, mt, poprank = 1,
    popunique = "MANIFESTMEANS", optimcontrol = list(tipredMissing = "sample")))
  slot <- sampled$ti_missing$parameter
  expect_equal(slot, length(x) + 1L)
  xs <- c(x, imputed$ti_imputed$mu + 0.3)
  at <- .pu_eval(sampled, xs)
  h <- 1e-5
  e <- replace(numeric(length(xs)), slot, h)
  fd <- (.pu_eval(sampled, xs + e)$value - .pu_eval(sampled, xs - e)$value) / (2 * h)
  expect_equal(as.numeric(at$gradient)[slot], fd, tolerance = 1e-6)
  expect_error(suppressWarnings(suppressMessages(ctFit(dt, mt, fit = FALSE,
    backend = "stan"))), "the stan backend does not", fixed = TRUE)
})

test_that("the gated floor on integrated deviations is the layer's, unit for unit", {
  skip_without_julia()
  # The gate reads the joint curvature's smallest eigenvalue, which the
  # profile's can exceed, and a unit it flags is scored in joint coordinates
  # with the filter at fixed deviations. The band is widened (hi 0.95) so that
  # units are flagged on this fixture at all, and the count is asserted so the
  # comparison cannot pass without the gated path running.
  m <- .pu_dyn_model()
  d <- .pu_dyn_data()
  old <- options(ctsem.laplace.hybrid = TRUE)
  on.exit(options(old), add = TRUE)
  spec <- .pu_spec(d, m, poprank = 1, popunique = "MANIFESTMEANS",
    optimcontrol = list(laplace_floor = "gated"))
  options(ctsem.laplace.hybrid = FALSE)
  layer <- .pu_spec(d, m, poprank = 1, popunique = "MANIFESTMEANS",
    optimcontrol = list(laplace_floor = "gated"))
  handle <- function(s) structure(s, class = c("ctJuliaModel", "ctFitModel"))
  objectives <- list(ctsem:::.ctJuliaObjective(handle(spec), hybrid = TRUE),
    ctsem:::.ctJuliaObjective(handle(layer), hybrid = TRUE))
  on.exit(ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_gated_adjoint!", TRUE),
    add = TRUE)
  previous <- lapply(objectives, function(o) ctsem:::.ctJuliaCall(
    "ContinuousTimeSEM.ctsem_set_laplace_floor!", o, "gated", hi = 0.95))
  on.exit(for (i in seq_along(objectives)) ctsem:::.ctJuliaCall(
    "ContinuousTimeSEM.ctsem_set_laplace_floor!", objectives[[i]], "gated",
    hi = previous[[i]]$hi), add = TRUE)
  # The count the last evaluation left, read as the field it is.
  gated <- function(o) as.integer(ctsem:::.ctJuliaCall("getfield", o,
    ctsem:::.ctJuliaEval(":gated_units")))
  # 10 and 18 units flagged at these two draws.
  for (seed in c(2, 4)) {
    set.seed(seed)
    x <- rnorm(length(.pu_names(spec)), 0, 1)
    full <- .pu_eval(layer, x)
    expect_gt(gated(objectives[[2]]), 0L)
    hybrid <- .pu_eval(spec, x)
    expect_equal(gated(objectives[[1]]), gated(objectives[[2]]))
    expect_equal(hybrid$value, full$value, tolerance = 1e-9)
    expect_equal(as.numeric(hybrid$gradient), as.numeric(full$gradient),
      tolerance = 1e-7)
    # The gated units' gradient by the reverse pass (the default) against
    # forward mode over the parameters, for the layer and the hybrid.
    for (s in list(layer, spec)) {
      ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_gated_adjoint!", FALSE)
      forward <- .pu_eval(s, x)
      ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_set_gated_adjoint!", TRUE)
      reverse <- .pu_eval(s, x)
      expect_equal(reverse$value, forward$value, tolerance = 1e-12)
      expect_equal(as.numeric(reverse$gradient), as.numeric(forward$gradient),
        tolerance = 1e-7)
    }
  }
})

test_that("quadrature integrates the filter's deviations at every node, and spans the shared effects alone", {
  skip_without_julia()
  # Under the hybrid a node's integrand is the members' likelihood with their
  # own deviations integrated exactly (the profile less logdet M / 2), so the
  # rule spans the one shared effect; on the layer it spans the effect and all
  # three deviations. Both converge on the same integral. The continuation on
  # the hybrid is that quadrature at its centre, its gradient is its value's,
  # and its Hessian is exact -- the own variances, which the filter reads,
  # seeded with the parameters it reads -- as central differences confirm.
  m <- .pu_dyn_model()
  d <- .pu_dyn_data()
  old <- options(ctsem.laplace.hybrid = TRUE)
  on.exit(options(old), add = TRUE)
  spec <- .pu_spec(d, m, poprank = 1, popunique = "MANIFESTMEANS")
  options(ctsem.laplace.hybrid = FALSE)
  layer <- .pu_spec(d, m, poprank = 1, popunique = "MANIFESTMEANS")
  handle <- function(s) structure(s, class = c("ctJuliaModel", "ctFitModel"))
  module <- ctsem:::.ctJuliaModule(spec$project)
  oh <- ctsem:::.ctJuliaObjective(handle(spec), hybrid = TRUE)
  ol <- ctsem:::.ctJuliaObjective(handle(layer), hybrid = TRUE)
  set.seed(5)
  x <- rnorm(length(.pu_names(spec)), 0, 0.3)
  xx <- ctsem:::.ctJuliaNumericVector(x)
  quad <- function(o, nodes) ctsem:::.ctBackendJuliaValue(
    module$ctsem_laplace_quadrature(o, xx, nodes = as.integer(nodes),
      contributions = TRUE))
  h9 <- quad(oh, 9)
  l9 <- quad(ol, 9)
  # Per subject, against the layer's 9^4-node rule.
  expect_equal(as.numeric(h9$subject_loglik), as.numeric(l9$subject_loglik),
    tolerance = 1e-6)
  expect_equal(quad(oh, 15)$value, h9$value, tolerance = 1e-10)

  cont <- module$ctsem_laplace_continuation(oh, xx, nodes = 5L, tolerance = -1)
  info <- ctsem:::.ctJuliaGet(module$ctsem_laplace_continuation_info(cont))
  expect_equal(as.integer(info$nflagged), 24L)
  expect_equal(as.integer(info$nwide), 0L)
  expect_equal(as.numeric(info$value), quad(oh, 5)$value, tolerance = 1e-10)
  value <- function(z) as.numeric(ctsem:::.ctJuliaGet(
    module$ctsem_laplace_continuation_evaluate(cont,
      ctsem:::.ctJuliaNumericVector(z), gradient = FALSE))$value)
  g <- as.numeric(ctsem:::.ctJuliaGet(module$ctsem_laplace_continuation_evaluate(
    cont, xx))$gradient)
  v <- rnorm(length(x))
  v <- v / sqrt(sum(v^2))
  h <- 1e-4
  expect_equal((value(x + h * v) - value(x - h * v)) / (2 * h), sum(g * v),
    tolerance = 1e-6)
  flagged <- ctsem:::.ctJuliaCall("ContinuousTimeSEM._continuation_flagged_hessian",
    cont, xx)
  expect_false(is.null(flagged))
  hess <- function(scheme) as.matrix(ctsem:::.ctBackendJuliaValue(
    module$ctsem_laplace_continuation_hessian(cont, xx, scheme = scheme)))
  expect_equal(hess("exact"), hess("central"), tolerance = 1e-5)
})

test_that("leave-one-out and linearised draws through the hybrid are the layer's", {
  skip_without_julia()
  # The deviations each filter integrates are drawn with the shared effect from
  # their joint Laplace posterior, in the layer's coordinates and order -- an
  # item a subject never answers from its prior -- so with one seed the draws,
  # their densities and every row's likelihood at them are the layer's. The
  # unit terms at parameter draws are the objective's, which is the layer's.
  m <- .pu_dyn_model()
  d <- .pu_dyn_data()
  old <- options(ctsem.laplace.hybrid = TRUE)
  on.exit(options(old), add = TRUE)
  spec <- .pu_spec(d, m, poprank = 1, popunique = "MANIFESTMEANS")
  options(ctsem.laplace.hybrid = FALSE)
  layer <- .pu_spec(d, m, poprank = 1, popunique = "MANIFESTMEANS")
  handle <- function(s) structure(s, class = c("ctJuliaModel", "ctFitModel"))
  module <- ctsem:::.ctJuliaModule(spec$project)
  objectives <- list(hybrid = ctsem:::.ctJuliaObjective(handle(spec), hybrid = TRUE),
    layer = ctsem:::.ctJuliaObjective(handle(layer), hybrid = TRUE))
  set.seed(9)
  x <- rnorm(length(.pu_names(spec)), 0, 0.3)
  draws <- lapply(objectives, function(o) ctsem:::.ctBackendJuliaValue(
    module$ctsem_laplace_effect_draws(o, ctsem:::.ctJuliaNumericVector(x), 12L,
      seed = 7L, scale = 1.2)))
  expect_true(all(as.logical(draws$hybrid$converged)))
  for (f in c("logq", "logprior", "llrow")) {
    expect_equal(as.numeric(draws$hybrid[[f]]), as.numeric(draws$layer[[f]]),
      tolerance = 1e-6, info = f)
  }
  at <- cbind(x, x + rnorm(length(x), 0, 0.05))
  terms <- lapply(objectives, function(o) ctsem:::.ctBackendJuliaValue(
    module$ctsem_unit_terms(o, ctsem:::.ctJuliaPut(at))))
  expect_equal(as.numeric(terms$hybrid$unit_loglik),
    as.numeric(terms$layer$unit_loglik), tolerance = 1e-8)
  # Subject parameters at parameter draws, each subject's modes linearised
  # around the estimate: the deviations' modes move with the parameters,
  # directly and through the shared effect's, as the layer's coordinates do.
  linear <- lapply(objectives, function(o) as.numeric(ctsem:::.ctBackendJuliaValue(
    module$ctsem_laplace_subject_values(o, ctsem:::.ctJuliaPut(t(at)),
      ctsem:::.ctJuliaNumericVector(x)))))
  expect_equal(linear$hybrid, linear$layer, tolerance = 1e-6)
})

test_that("reporting through the hybrid is the layer's: population, subjects, modes, Kalman", {
  skip_without_julia()
  # Each subject's parameters carry its filter-integrated deviations at their
  # conditional modes, the population sds their own variances, the modes
  # report their joint posterior sd, and the Kalman output is a plug-in filter
  # at those parameters: all what the layer holding the deviations as
  # coordinates reports, on a model with shared effects on the dynamics and on
  # the diagonal (rank 0) one.
  old <- options(ctsem.laplace.hybrid = TRUE)
  on.exit(options(old), add = TRUE)
  d <- .pu_dyn_data()
  handle <- function(s) structure(s, class = c("ctJuliaModel", "ctFitModel"))
  cases <- list(list(m = .pu_dyn_model(), rank = 1, pu = "MANIFESTMEANS"),
    list(m = .pu_model(), rank = 0, pu = TRUE))
  set.seed(21)
  for (case in cases) {
    options(ctsem.laplace.hybrid = TRUE)
    hy <- .pu_spec(d, case$m, poprank = case$rank, popunique = case$pu)
    options(ctsem.laplace.hybrid = FALSE)
    la <- .pu_spec(d, case$m, poprank = case$rank, popunique = case$pu)
    expect_false(is.null(hy$laplace$hybrid))
    x <- rnorm(length(.pu_names(hy)), 0, 0.3)
    module <- ctsem:::.ctJuliaModule(hy$project)
    report <- function(s) {
      o <- ctsem:::.ctJuliaObjective(handle(s), hybrid = TRUE)
      xx <- ctsem:::.ctJuliaNumericVector(x)
      list(
        population = ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_population(o,
          ctsem:::.ctJuliaPut(matrix(x, 1)), 1L)),
        subjects = as.numeric(ctsem:::.ctBackendJuliaValue(
          module$ctsem_laplace_subject_values(o, xx))),
        modes = ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_modes(o, xx, 1L)),
        kalman = ctsem:::.ctBackendJuliaValue(module$ctsem_kalman(o, xx,
          subject_matrices = FALSE)))
    }
    rh <- report(hy)
    rl <- report(la)
    expect_equal(as.numeric(rh$population$sd), as.numeric(rl$population$sd),
      tolerance = 1e-10)
    expect_equal(as.numeric(rh$population$correlation),
      as.numeric(rl$population$correlation), tolerance = 1e-10)
    expect_equal(rh$subjects, rl$subjects, tolerance = 1e-7)
    expect_equal(as.numeric(rh$modes$raw), as.numeric(rl$modes$raw), tolerance = 1e-7)
    expect_equal(as.numeric(rh$modes$raw_sd), as.numeric(rl$modes$raw_sd),
      tolerance = 1e-7)
    for (f in c("eta", "etacov", "y", "ycov", "llrow")) {
      expect_equal(as.numeric(rh$kalman[[f]]), as.numeric(rl$kalman[[f]]),
        tolerance = 1e-7, info = f)
    }
  }
})

test_that("sequential predictions are the fitted model's own, row by row", {
  skip_without_julia()
  # Each row's prediction from the subject's earlier rows alone, the effects at
  # their Laplace posterior given those rows. Checked: the rows' log predictive
  # densities sum to the objective (no priors, so it is the subjects' terms
  # alone); the last row's updated values, conditional on the whole record,
  # are the plug-in filter's at the fitted modes; a subject's first prediction,
  # with nothing observed, adds the population covariance of the intercepts to
  # the filter's own; and the hybrid and the layer agree.
  old <- options(ctsem.laplace.hybrid = TRUE)
  on.exit(options(old), add = TRUE)
  d <- .pu_dyn_data()
  handle <- function(s) structure(s, class = c("ctJuliaModel", "ctFitModel"))
  set.seed(31)
  cases <- list(list(m = .pu_model(), linear = TRUE),
    list(m = .pu_dyn_model(), linear = FALSE))
  for (case in cases) {
    out <- list()
    x <- NULL
    for (h in c(TRUE, FALSE)) {
      options(ctsem.laplace.hybrid = h)
      spec <- .pu_spec(d, case$m, poprank = 1, popunique = "MANIFESTMEANS",
        priors = FALSE)
      if (is.null(x)) x <- rnorm(length(.pu_names(spec)), 0, 0.3)
      module <- ctsem:::.ctJuliaModule(spec$project)
      o <- ctsem:::.ctJuliaObjective(handle(spec), hybrid = TRUE)
      xx <- ctsem:::.ctJuliaNumericVector(x)
      seq <- ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_sequential(o, xx))
      value <- ctJuliaEvaluate(handle(spec), x, gradient = FALSE, cores = 1L)$value
      expect_equal(sum(seq$llrow), value, tolerance = 1e-8)
      kal <- ctsem:::.ctBackendJuliaValue(module$ctsem_kalman(o, xx,
        subject_matrices = FALSE))
      subject <- as.integer(seq$subject)
      last <- as.integer(tapply(seq_along(subject), subject, max))
      expect_equal(seq$eta[2, last, , drop = FALSE], kal$eta[2, last, , drop = FALSE],
        tolerance = 1e-7)
      expect_equal(seq$y[2, last, , drop = FALSE], kal$y[2, last, , drop = FALSE],
        tolerance = 1e-7)
      if (case$linear && h) {
        first <- as.integer(tapply(seq_along(subject), subject, min))
        base <- ctsem:::.ctBackendJuliaValue(module$ctsem_kalman(o, xx,
          subject_matrices = FALSE, from_level = 2L))
        popcov <- as.matrix(ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_popcov(o,
          xx, 1L)))
        slope <- spec$laplace$hybrid$slope
        rows <- spec$laplace$hybrid$manifest
        added <- matrix(0, 3, 3)
        added[rows, rows] <- diag(slope, 3) %*% popcov %*% diag(slope, 3)
        for (r in first) {
          expect_equal(seq$ycov[1, r, , ], base$ycov[1, r, , ] + added, tolerance = 1e-8)
        }
      }
      out[[as.character(h)]] <- seq
    }
    for (f in c("eta", "etacov", "y", "ycov", "llrow")) {
      expect_equal(out[["TRUE"]][[f]], out[["FALSE"]][[f]], tolerance = 1e-7, info = f)
    }
  }
})

test_that("a sequential prediction does not see the subject's later data", {
  skip_without_julia()
  # Subject 1's observations after time 2 withheld: its predicted and updated
  # values up to then are unchanged. The plug-in filter, at modes from the
  # whole record, moves -- which is what makes this a test of anything.
  old <- options(ctsem.laplace.hybrid = TRUE)
  on.exit(options(old), add = TRUE)
  handle <- function(s) structure(s, class = c("ctJuliaModel", "ctFitModel"))
  d <- .pu_dyn_data()
  cut <- d
  cut[cut$id == 1 & cut$time > 2, c("Y1", "Y2", "Y3")] <- NA
  m <- .pu_dyn_model()
  set.seed(41)
  x <- NULL
  report <- function(dd) {
    spec <- .pu_spec(dd, m, poprank = 1, popunique = "MANIFESTMEANS")
    if (is.null(x)) x <<- rnorm(length(.pu_names(spec)), 0, 0.3)
    module <- ctsem:::.ctJuliaModule(spec$project)
    o <- ctsem:::.ctJuliaObjective(handle(spec), hybrid = TRUE)
    xx <- ctsem:::.ctJuliaNumericVector(x)
    list(seq = ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_sequential(o, xx)),
      plug = ctsem:::.ctBackendJuliaValue(module$ctsem_kalman(o, xx,
        subject_matrices = FALSE)))
  }
  full <- report(d)
  part <- report(cut)
  rows <- which(full$seq$subject == 1)[d$time[d$id == 1] <= 2]
  at <- function(a) if (length(dim(a)) == 3L) a[, rows, , drop = FALSE] else
    a[, rows, , , drop = FALSE]
  for (f in c("eta", "etacov", "y", "ycov")) {
    expect_equal(at(full$seq[[f]]), at(part$seq[[f]]), tolerance = 1e-9, info = f)
  }
  expect_false(isTRUE(all.equal(full$plug$y[1, rows, , drop = FALSE],
    part$plug$y[1, rows, , drop = FALSE], tolerance = 1e-9)))
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
  # The fit integrates the own variances in its filter, and its reports read
  # them from that objective: the same as the layer holding them as
  # coordinates reports, which the fit with its hybrid record removed uses.
  expect_false(is.null(fit$model_spec$laplace$hybrid))
  layer <- fit
  layer$model_spec$laplace$hybrid <- NULL
  expect_equal(summary(layer)$popsd, s$popsd, tolerance = 1e-8)
  expect_equal(ctSubjectPars(layer), ctSubjectPars(fit), tolerance = 1e-6)
  quiet <- function(x) suppressWarnings(suppressMessages(x))
  kh <- quiet(ctKalman(fit, subjects = 1:3))
  kl <- quiet(ctKalman(layer, subjects = 1:3))
  numeric <- vapply(kh, is.numeric, logical(1))
  expect_equal(kh[numeric], kl[numeric], tolerance = 1e-6)
})
