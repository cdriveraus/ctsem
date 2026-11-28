# ctFit(intoverpop = 'variational'): each unit's random effects integrated by
# a lower bound fitted at every evaluation (variational.jl in the engine,
# `.ctBackendVariational()`).
#
# What is pinned here: on a model where every unit's integrand is Gaussian the
# bound is the integral, so the fit is Laplace's (exact there); on a
# nonlinear one the fit reports the random-effects check's integral with the
# bound below it and the check's step from the estimate, its estimates are
# those of the refined fit, which integrates every unit numerically, and
# without the check it reports the bound; and the route is refused where it
# cannot apply.

.vf_linear_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), MANIFESTMEANS = matrix("mmean"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE
  model
}

.vf_linear_data <- function(nsubjects = 25, nobs = 6) {
  set.seed(20260827)
  drift <- -0.4; diffusion <- 0.6
  do.call(rbind, lapply(seq_len(nsubjects), function(i) {
    intercept <- stats::rnorm(1, 1.5, 0.9)
    state <- stats::rnorm(1, 0, 0.5)
    out <- numeric(nobs)
    for (t in seq_len(nobs)) {
      if (t > 1) {
        decay <- exp(drift)
        state <- decay * state +
          stats::rnorm(1, 0, sqrt(diffusion^2 / (-2 * drift) * (1 - decay^2)))
      }
      out[t] <- state + intercept + stats::rnorm(1, 0, 0.3)
    }
    data.frame(id = i, time = seq_len(nobs) - 1, Y1 = out)
  }))
}

# One random effect, on a `-log1p_exp` DRIFT, simulated here rather than by
# ctGenerate, whose draw stream moves.
.vf_nonlinear_model <- function() {
  model <- suppressMessages(ctModel(silent = TRUE, type = "ct", CINT = 0,
    MANIFESTMEANS = 0, LAMBDA = matrix(1), T0MEANS = matrix(0),
    DRIFT = "drift|-log1p_exp(-param)|TRUE"))
  model$pars$indvarying <- model$pars$param %in% "drift"
  model
}

.vf_nonlinear_data <- function(seed = 3L, nsubjects = 40L, ntimes = 8L) {
  set.seed(seed)
  drift <- -log1p(exp(-stats::rnorm(nsubjects, 1, 1)))
  rows <- lapply(seq_len(nsubjects), function(i) {
    a <- drift[i]; decay <- exp(a)
    innovation <- sqrt(0.25 * (exp(2 * a) - 1) / (2 * a))
    latent <- numeric(ntimes); latent[1] <- stats::rnorm(1, 0, 1)
    for (t in seq_len(ntimes - 1L)) latent[t + 1L] <- decay * latent[t] +
      stats::rnorm(1, 0, innovation)
    data.frame(id = i, time = seq_len(ntimes) - 1L,
      Y1 = latent + stats::rnorm(ntimes, 0, 0.3))
  })
  do.call(rbind, rows)
}

.vf_cache <- new.env(parent = emptyenv())
.vf_fit <- function(which, route) {
  key <- paste(which, route)
  if (!exists(key, envir = .vf_cache, inherits = FALSE)) {
    data <- if (identical(which, "linear")) .vf_linear_data() else .vf_nonlinear_data()
    model <- if (identical(which, "linear")) .vf_linear_model() else .vf_nonlinear_model()
    set.seed(1)
    assign(key, suppressWarnings(suppressMessages(ctFit(data, model,
      backend = "julia", intoverpop = route, cores = 2))), envir = .vf_cache)
  }
  get(key, envir = .vf_cache)
}

test_that("on Gaussian units the variational fit is Laplace's", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_fit("linear", "variational")
  lf <- .vf_fit("linear", "laplace")
  expect_identical(vf$args$resolved$intoverpop, "variational")
  expect_identical(vf$estimate$loglik_method, "variational")
  expect_equal(vf$estimate$raw, lf$estimate$raw, tolerance = 1e-3)
  expect_equal(vf$estimate$loglik, lf$estimate$loglik, tolerance = 1e-3)
  expect_equal(vf$estimate$loglik_lowerbound, vf$estimate$loglik, tolerance = 1e-3)
  expect_equal(vf$estimate$se, lf$estimate$se, tolerance = 0.02)
})

test_that("on a nonlinear model the bound lies below the reported log likelihood", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_fit("nonlinear", "variational")
  qf <- .vf_fit("nonlinear", "quadrature")
  e <- vf$estimate
  expect_identical(e$loglik_method, "variational")
  expect_true(is.finite(e$loglik) && is.finite(e$loglik_lowerbound))
  expect_lte(e$loglik_lowerbound, e$loglik + 4 * e$loglik_se)
  expect_equal(ctsem:::.ctFitOptimisedLoglik(vf), e$loglik_lowerbound)
  # The refined fit integrates every flagged unit numerically; the variational
  # estimate lies within a fifth of a standard error of it.
  expect_lt(max(abs(vf$estimate$raw - qf$estimate$raw) / qf$estimate$se), 0.2)
  expect_equal(e$loglik, qf$estimate$loglik, tolerance = 0.1)
  # The check's step from the estimate towards the integral's maximum, in the
  # bound's standard errors: small where the estimate is the refined fit's.
  expect_true(vf$variational$check$status %in% c("passed", "unresolved", "flagged"))
  expect_lt(vf$variational$check$max_shift, 0.2)
})

test_that("a variational fit's interval check and 'checked' run", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_fit("nonlinear", "variational")
  check <- vf$uncertainty$details$normalcheck
  expect_null(check$failed)
  expect_gt(check$draws, 0)
  # 'checked' takes no coreset of an objective the engine cannot take one of.
  u <- suppressWarnings(suppressMessages(ctFitUncertainty(vf, "checked", cores = 2)))
  expect_true(all(is.finite(u$estimate$se)))
  # The joint sample's target is the Laplace objective of the same
  # specification; the bound's marginal is no sampler's.
  lf <- vf
  lf$model_spec$laplace$variational <- NULL
  joint <- ctsem:::.ctBackendSampleObjective(vf, list(marginal = FALSE))
  x <- as.numeric(vf$estimate$raw)
  module <- ctsem:::.ctJuliaModule(vf$model_spec$project)
  value <- function(o) as.numeric(ctsem:::.ctJuliaGet(module$ctsem_evaluate(o,
    ctsem:::.ctJuliaNumericVector(x), gradient = FALSE))$value)
  expect_equal(value(joint), value(ctsem:::.ctJuliaObjective(lf)))
  expect_error(ctsem:::.ctBackendSampleObjective(vf, list(marginal = TRUE)),
    "marginal cannot be sampled")
})

# The interval check on a variational fit solves its draws to a rough
# tolerance (`.ctVariationalCheckTolerance`), and the fit's own is back once
# it has run: what the objective's tolerance is while the check runs, read by
# a stand-in for the check, and after.
test_that("a variational fit's interval check solves its draws roughly", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_fit("nonlinear", "variational")
  objective <- ctsem:::.ctJuliaObjective(vf)
  tolerance <- function() {
    old <- as.numeric(ctsem:::.ctJuliaCall(
      "ContinuousTimeSEM.ctsem_variational_tolerance!", objective, 1))
    ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_variational_tolerance!", objective, old)
    old
  }
  expect_equal(tolerance(), 1e-8)
  during <- NULL
  testthat::local_mocked_bindings(.ctOptimNormalCheck = function(gradients, est, cov,
    ...) {
    during <<- tolerance()
    list(draws = 0L, unusable = 0L, directions = matrix(0, length(est), 0L),
      slope = numeric(), shape = numeric(), widths = numeric(),
      parameters = integer(), centre = as.numeric(est))
  })
  suppressWarnings(suppressMessages(ctFitUncertainty(vf, cores = 2)))
  expect_equal(during, ctsem:::.ctVariationalCheckTolerance)
  expect_equal(tolerance(), 1e-8)
})

test_that("a variational fit is refined by quadrature from its own estimate", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_fit("nonlinear", "variational")
  qf <- .vf_fit("nonlinear", "quadrature")
  rf <- suppressWarnings(suppressMessages(ctFitUncertainty(vf,
    intoverpop = "quadrature", cores = 2)))
  corr <- rf$laplace$correction
  expect_identical(corr$status, "continued")
  expect_true(corr$applied)
  # The climb started at the variational estimate, on the Laplace objective
  # with the units the check names against Laplace there flagged.
  expect_equal(corr$laplace_estimate, vf$estimate$raw)
  expect_null(rf$model_spec$laplace$variational)
  expect_gt(length(corr$check$carriers), 0)
  expect_true(all(corr$check$carriers %in% corr$flagged_units))
  expect_identical(rf$estimate$loglik_method, "quadrature")
  # It ends where ctFit(intoverpop = 'quadrature') ends from Laplace's
  # optimum, up to the units the screen flags at the two starts.
  expect_lt(max(abs(rf$estimate$raw - qf$estimate$raw) / qf$estimate$se), 0.1)
  expect_lt(abs(rf$estimate$loglik - qf$estimate$loglik), 0.1)
})

test_that("without the random-effects check a variational fit reports its bound", {
  skip_on_cran()
  skip_without_julia()
  set.seed(1)
  vf <- suppressWarnings(suppressMessages(ctFit(.vf_nonlinear_data(),
    .vf_nonlinear_model(), backend = "julia", intoverpop = "variational", cores = 2,
    optimcontrol = list(checks = FALSE))))
  e <- vf$estimate
  expect_identical(e$loglik_method, "variational_bound")
  expect_identical(e$loglik, e$loglik_lowerbound)
  expect_null(vf$variational)
  expect_equal(e$loglik, .vf_fit("nonlinear", "variational")$estimate$loglik_lowerbound,
    tolerance = 1e-4)
})

test_that("leave-one-out and the score tests read a variational fit's bound", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_fit("nonlinear", "variational")
  est <- as.numeric(vf$estimate$raw)
  # The units' terms are their bounds, and a subject with every row withheld
  # scores nothing.
  ut <- ctsem:::.ctBackendUnitTerms(vf$model_spec, est)
  expect_equal(sum(ut$unit_loglik), vf$estimate$loglik_lowerbound, tolerance = 1e-6)
  training <- vf$model_spec$data
  training[training$id == unique(training$id)[3], "Y1"] <- NA
  tf <- ctsem:::.ctFitReplaceData(vf, training)
  emptied <- ctsem:::.ctBackendUnitTerms(tf$model_spec, est, subjects = 3L)
  expect_equal(emptied$unit_loglik[3, 1], 0, tolerance = 1e-8)
  l <- suppressWarnings(suppressMessages(ctLOO(vf, folds = 4, cores = 2,
    casewiseApproximation = TRUE)))
  expect_true(all(is.finite(l$outsampleLogLikSubject)))
  expect_match(l$scoring, "variational bound")
  # A subject's row scores sum to its score on the bound, up to the inner
  # solves' tolerance.
  R <- ctsem:::.ctBackendRowScores(vf, est, diag(length(est)), 2)
  S <- ctsem:::.ctBackendScoreMatrix(vf, est, prior = FALSE)
  expect_lt(max(abs(rowsum(R, ctsem:::.ctFitRowSubject(vf)) - S)), 1e-3)
  # A sample takes the joint; ctLaplaceCheck reads the integral against the
  # bound.
  expect_identical(ctsem:::.ctBackendSampleRoute("laplace", requested = "variational"),
    "none")
  chk <- suppressWarnings(suppressMessages(ctLaplaceCheck(vf, cores = 2)))
  expect_identical(chk$approximation, "variational")
  expect_equal(chk$laplace, vf$estimate$loglik_lowerbound, tolerance = 1e-6)
  expect_equal(chk$gap, chk$quadrature - chk$laplace)
})

# Persons in studies: a manifest mean that varies at both levels (a Gaussian
# integrand, where the bound is the integral), and with `drift` a person-level
# drift through `-log1p_exp` beside it (not Gaussian).
.vf_nested_model <- function(drift = FALSE) {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", silent = TRUE,
    id = c("id", "study"), manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), DRIFT = if (drift) "drift|-log1p_exp(-param)|TRUE" else
      matrix(-0.5), DIFFUSION = matrix(0.5), MANIFESTMEANS = matrix("mm"),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0))))
  m$pars$indvarying <- m$pars$param %in% "mm" | (drift & m$pars$param %in% "drift")
  if (drift) m$pars$indvarying[m$pars$param %in% "drift"] <- "id"
  m
}

.vf_nested_data <- function(npersons = 40L, nstudies = 8L, ntimes = 8L, seed = 7L) {
  set.seed(seed)
  study <- rep(seq_len(nstudies), length.out = npersons)
  smean <- stats::rnorm(nstudies, 1, 0.7)
  rows <- lapply(seq_len(npersons), function(i) {
    a <- -log1p(exp(-stats::rnorm(1, 0.7, 0.6)))
    decay <- exp(a)
    innovation <- sqrt(0.25 * (exp(2 * a) - 1) / (2 * a))
    latent <- numeric(ntimes)
    latent[1] <- stats::rnorm(1, 0, 1)
    for (t in seq_len(ntimes - 1L)) latent[t + 1L] <- decay * latent[t] +
      stats::rnorm(1, 0, innovation)
    data.frame(id = i, study = study[i], time = seq_len(ntimes) - 1L,
      Y1 = smean[study[i]] + stats::rnorm(1, 0, 0.5) + latent +
        stats::rnorm(ntimes, 0, 0.3))
  })
  do.call(rbind, rows)
}

.vf_nested_fit <- function(drift, route) {
  key <- paste("nested", drift, route)
  if (!exists(key, envir = .vf_cache, inherits = FALSE)) {
    set.seed(1)
    assign(key, suppressWarnings(suppressMessages(ctFit(.vf_nested_data(),
      .vf_nested_model(drift), backend = "julia", intoverpop = route, cores = 2))),
      envir = .vf_cache)
  }
  get(key, envir = .vf_cache)
}

test_that("persons in studies: the tree family is Laplace's on a Gaussian integrand", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_nested_fit(FALSE, "variational")
  lf <- .vf_nested_fit(FALSE, "laplace")
  expect_identical(vf$estimate$loglik_method, "variational")
  expect_equal(vf$estimate$raw, lf$estimate$raw, tolerance = 1e-3)
  expect_equal(vf$estimate$loglik_lowerbound, vf$estimate$loglik, tolerance = 1e-3)
})

test_that("persons in studies with a random drift: the bound lies below the integral", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_nested_fit(TRUE, "variational")
  e <- vf$estimate
  expect_identical(e$loglik_method, "variational")
  expect_true(is.finite(e$loglik) && is.finite(e$loglik_lowerbound))
  expect_lte(e$loglik_lowerbound, e$loglik + 4 * e$loglik_se + e$loglik_bound)
})

test_that("a variational fit's subject effects are its approximations' means", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_fit("nonlinear", "variational")
  est <- as.numeric(vf$estimate$raw)
  V <- ctsem:::.ctJuliaObjective(vf)
  # Each person's posterior mean of their drift effect by a dense grid of
  # their exact integrand, against the mean of their approximation and the
  # mode Laplace reports, in posterior sds.
  ctsem:::.ctJuliaEval(paste(
    "function _vf_test_means(V, x)",
    "  C = ContinuousTimeSEM",
    "  lap = V.laplace",
    "  lam = C._vf_solved(V, x)",
    "  theta = C._laplace_extend(lap, collect(Float64, x))",
    "  Ls = [Matrix{Float64}(M) for M in C._laplace_popchols(theta, lap.spec)]",
    "  aws = C._laplace_workspace!(lap, Float64, length(theta))",
    "  grid = collect(-8.0:0.005:8.0)",
    "  out = zeros(length(lam), 4)",
    "  for U in eachindex(lam)",
    "    g = [C._laplace_unit_objective_gradient(lap, U, theta, Ls, [u], aws).value for u in grid]",
    "    w = exp.(g .- maximum(g)); w ./= sum(w)",
    "    ex = sum(w .* grid)",
    "    out[U, :] = [C._vf_unit_mean(V.layouts[U], lam[U])[1],",
    "      C._laplace_newton_unit_mode(lap, U, theta, Ls, aws, zeros(1)).u[1], ex,",
    "      sqrt(sum(w .* (grid .- ex) .^ 2))]",
    "  end",
    "  out",
    "end", sep = "\n"))
  r <- matrix(as.numeric(ctsem:::.ctJuliaCall("Main._vf_test_means", V,
    ctsem:::.ctJuliaNumericVector(est))), ncol = 4)
  expect_lt(max(abs(r[, 1] - r[, 3]) / r[, 4]), 0.01)
  expect_gt(stats::median(abs(r[, 2] - r[, 3]) / r[, 4]), 0.05)
  # Each subject's raw drift is the population's shifted by that mean.
  module <- ctsem:::.ctJuliaModule(vf$model_spec$project)
  sv <- as.matrix(ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_subject_values(V,
    ctsem:::.ctJuliaNumericVector(est), ti_effects = FALSE)))
  drift <- match("drift", ctsem:::.ctBackendRawParameterNames(vf, length(est)))
  expect_equal(cor(sv[, drift], r[, 1]), 1, tolerance = 1e-6)
  # And the readers of the effects run on them.
  expect_equal(dim(ctSubjectPars(vf)), c(1L, 40L, 1L))
  expect_true(nrow(suppressMessages(ctKalman(vf, subjects = 1:2))) > 0)
  ex <- ctExtract(vf, subjectMatrices = TRUE, nsamples = 10)
  expect_false(is.null(ex))
})

# A person-specific loading whose sign the data barely say: each person's
# posterior has two maxima, one each side of zero (test-julia-laplace-modes.R's
# sign fixture), so each unit's bound has an optimum near each.
.vf_sign_data <- function() {
  set.seed(3)
  n <- 12; nt <- 8
  d <- data.frame(id = rep(1:n, each = nt), time = rep(0:(nt - 1), n))
  lambda <- 1.5 * sample(c(-1, 1), n, TRUE)
  eta <- unlist(lapply(1:n, function(i) {
    e <- numeric(nt)
    e[1] <- stats::rnorm(1)
    for (t in 2:nt) e[t] <- 0.6 * e[t - 1] + stats::rnorm(1, 0, 0.8)
    e
  }))
  d$Y1 <- lambda[d$id] * eta + stats::rnorm(n * nt, 0, 0.3)
  d
}

.vf_sign_model <- function() {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", id = "id",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix("lambda"),
    DRIFT = matrix(-0.5), DIFFUSION = matrix(0.8), MANIFESTVAR = matrix(0.3),
    MANIFESTMEANS = matrix(0), T0MEANS = matrix(0), T0VAR = matrix(1),
    CINT = matrix(0))))
  m$pars$indvarying <- m$pars$param %in% "lambda"
  m$pars$transform[m$pars$param %in% "lambda"] <- "param"
  m
}

test_that("a variational fit keeps its approximations, and a stage that reads it holds them", {
  skip_on_cran()
  skip_without_julia()
  set.seed(1)
  vf <- suppressWarnings(suppressMessages(ctFit(.vf_sign_data(), .vf_sign_model(),
    backend = "julia", intoverpop = "variational", cores = 2)))
  e <- vf$estimate
  modes <- ctsem:::.ctJuliaFitModes(vf)
  expect_identical(modes$units, 1:12)
  spec <- structure(vf$model_spec, class = c("ctJuliaModel", "ctFitModel"))
  V <- ctsem:::.ctJuliaObjective(spec)
  value <- function(x) as.numeric(ctJuliaEvaluate(spec, x, gradient = FALSE,
    cores = 1L)$value)
  at <- value(e$raw)
  # Outside a run an evaluation keeps nothing, so a line through the estimate
  # reads the same whichever way it is walked, and without a jump.
  grid <- e$raw[1] + seq(-2, 2, by = 0.25)
  out <- vapply(grid, function(a) value(c(a, e$raw[2])), numeric(1))
  back <- rev(vapply(rev(grid), function(a) value(c(a, e$raw[2])), numeric(1)))
  expect_identical(out, back)
  bend <- abs(diff(out, differences = 2))
  expect_lt(max(bend), 2 * stats::median(bend))
  # An objective that has forgotten them (as a new session's has), read first
  # at the mirrored loading, follows other optima; the fit's held
  # approximations give the fit's value back.
  reset <- function() invisible(ctsem:::.ctJuliaCall(
    "ContinuousTimeSEM.ctsem_variational_reset!", V))
  reset()
  invisible(value(c(-e$raw[1], e$raw[2])))
  held <- local({
    expect_identical(ctsem:::.ctJuliaHoldModes(spec, modes), 12L)
    value(e$raw)
  })
  expect_equal(held, at, tolerance = 1e-6)
  # An objective built for the stored fit, as a new session builds it, starts
  # from the fit's approximations, so read first at the mirrored loading it
  # still gives the fit's value at the estimate.
  cache <- get(".ct_julia_cache", envir = asNamespace("ctsem"))
  cache$objectives <- new.env(parent = emptyenv())
  invisible(ctsem:::.ctJuliaObjective(vf))
  invisible(value(c(-e$raw[1], e$raw[2])))
  expect_equal(value(e$raw), at, tolerance = 1e-6)
  # The reported integral is at least the bound wherever the draws, centred on
  # one of Laplace's maxima, would have missed the other's mass.
  expect_gte(e$loglik, e$loglik_lowerbound - 3 * e$loglik_se)
})

test_that("intoverpop = 'variational' is refused where it cannot apply", {
  skip_on_cran()
  expect_error(suppressWarnings(suppressMessages(ctFit(.vf_linear_data(5, 3),
    .vf_linear_model(), backend = "stan", intoverpop = "variational"))),
    "needs backend = 'julia'")
  expect_error(suppressWarnings(suppressMessages(ctFit(.vf_linear_data(5, 3),
    .vf_linear_model(), backend = "julia", intoverpop = "variational",
    optimize = FALSE))), "fits by optimizing")
})
