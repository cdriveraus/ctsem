# The quadrature correction every Laplace fit gets by default
# (`optimcontrol$laplace_correct`; `.ctLaplaceAutoCorrect()` in
# R/ctBackendLaplaceCorrect.R, `ctsem_laplace_autocorrect` in the engine).
#
# What is pinned here: that a model on which Laplace is exact is left alone to
# the bit, that a nonlinear one is moved towards the quadrature optimum and
# says so, that switching it off is the old fit, that the post-hoc functions do
# not apply it twice, that each request it cannot honour is refused by name,
# and that a fit without the record -- saved before it existed -- reads as
# uncorrected.

# A random effect on an identity-transformed MANIFESTMEANS: Laplace is exact.
.ac_linear_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), MANIFESTMEANS = matrix("mmean"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE
  model
}

.ac_linear_data <- function(nsubjects = 25, nobs = 6) {
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

# One random effect, on a `-log1p_exp` DRIFT: Laplace is not exact, and with
# one effect per unit the 5-node quadrature is 5 filter passes a unit.
# Simulated here rather than by ctGenerate, whose draw stream moves.
.ac_nonlinear_model <- function() {
  model <- suppressMessages(ctModel(silent = TRUE, type = "ct", CINT = 0,
    MANIFESTMEANS = 0, LAMBDA = matrix(1), T0MEANS = matrix(0),
    DRIFT = "drift|-log1p_exp(-param)|TRUE"))
  model$pars$indvarying <- model$pars$param %in% "drift"
  model
}

.ac_nonlinear_data <- function(seed = 3L, nsubjects = 40L, ntimes = 8L) {
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

.ac_cache <- new.env(parent = emptyenv())
# Each pair shares a seed, so everything before the correction -- start,
# optimiser, Hessian, the draws -- is the same computation in both.
.ac_fits <- function(which) {
  key <- paste0(which, "_fits")
  if (!exists(key, envir = .ac_cache, inherits = FALSE)) {
    data <- if (identical(which, "linear")) .ac_linear_data() else .ac_nonlinear_data()
    model <- if (identical(which, "linear")) .ac_linear_model() else .ac_nonlinear_model()
    fitwith <- function(...) {
      set.seed(11)
      suppressMessages(ctFit(data, model, backend = "julia",
        intoverpop = "laplace", cores = 1,
        optimcontrol = list(finishsamples = 100, ...)))
    }
    assign(key, list(on = fitwith(), off = fitwith(laplace_correct = FALSE)),
      envir = .ac_cache)
  }
  get(key, envir = .ac_cache, inherits = FALSE)
}

.ac_quadrature <- function(fit, at) {
  module <- .ctJuliaModule(fit$model_spec$project)
  JuliaConnectoR::juliaGet(module$ctsem_laplace_quadrature(.ctJuliaObjective(fit),
    .ctJuliaNumericVector(as.numeric(at)), nodes = 5L))$value
}

test_that("a model on which Laplace is exact is screened and left alone", {
  skip_without_julia()
  fits <- .ac_fits("linear")
  corr <- fits$on$laplace$correction
  expect_identical(corr$status, "exact")
  expect_false(corr$applied)
  expect_lt(corr$screen, corr$tolerance)
  # Screened, so no step was computed at all.
  expect_identical(corr$steps, 0L)
  expect_equal(as.numeric(fits$on$estimate$raw), as.numeric(fits$off$estimate$raw),
    tolerance = 1e-10)
  expect_equal(fits$on$estimate$loglik, fits$off$estimate$loglik, tolerance = 1e-10)
  expect_null(fits$on$estimate$loglik_method)
  expect_identical(fits$off$laplace$correction$status, "off")
  expect_false(any(grepl("corrected by quadrature", capture.output(print(fits$on)))))
})

test_that("a nonlinear random effect is corrected towards the quadrature optimum", {
  skip_without_julia()
  fits <- .ac_fits("nonlinear")
  on <- fits$on; off <- fits$off
  corr <- on$laplace$correction
  expect_identical(corr$status, "corrected")
  expect_true(corr$applied)
  expect_gt(corr$screen, corr$tolerance)
  expect_gte(corr$steps, 1L)
  # Where it started is the uncorrected fit, exactly.
  expect_equal(corr$laplace_estimate, as.numeric(off$estimate$raw), tolerance = 1e-10)
  expect_equal(corr$loglik_laplace, off$estimate$loglik, tolerance = 1e-10)
  expect_equal(as.numeric(on$estimate$raw), corr$laplace_estimate +
    unname(corr$delta), tolerance = 1e-12)
  expect_gt(max(abs(corr$delta)), 0)

  # The quadrature objective is higher at the corrected point, and the
  # corrected point is nearer the quadrature optimum than the Laplace one.
  expect_gt(.ac_quadrature(on, on$estimate$raw), .ac_quadrature(on, off$estimate$raw))
  module <- .ctJuliaModule(on$model_spec$project)
  refined <- as.numeric(JuliaConnectoR::juliaGet(module$ctsem_laplace_refine(
    .ctJuliaObjective(on), .ctJuliaNumericVector(corr$laplace_estimate),
    nodes = 5L, maxiter = 100L))$minimizer)
  expect_lt(sqrt(sum((on$estimate$raw - refined)^2)),
    sqrt(sum((off$estimate$raw - refined)^2)))

  # What the fit reports: the quadrature log likelihood at the corrected point,
  # with the Laplace one kept beside it, and the per-subject terms summing to it.
  expect_identical(on$estimate$loglik_method, "quadrature")
  expect_equal(on$estimate$loglik, corr$loglik_quadrature)
  expect_equal(on$estimate$loglik_laplace, off$estimate$loglik, tolerance = 1e-10)
  expect_equal(sum(on$estimate$subject_loglik), on$estimate$loglik, tolerance = 1e-8)
  expect_equal(on$estimate$logposterior, .ac_quadrature(on, on$estimate$raw),
    tolerance = 1e-8)

  # Recentred, not reshaped: the covariance is the Laplace curvature, still
  # said to be at the Laplace optimum, and the draws are the uncorrected fit's
  # moved by the step.
  expect_equal(as.matrix(on$estimate$cov), as.matrix(off$estimate$cov))
  expect_equal(as.numeric(on$uncertainty$evaluated_at), corr$laplace_estimate)
  shift <- on$estimate$rawposterior - off$estimate$rawposterior
  expect_equal(unname(shift), matrix(corr$delta, nrow(shift), ncol(shift),
    byrow = TRUE), tolerance = 1e-10)
  expect_equal(length(corr$delta_se), length(on$estimate$raw))
  expect_true(is.integer(corr$dropped_directions))
  expect_true(all(is.finite(corr$seconds)))
})

test_that("print says so only when the estimate moved materially", {
  skip_without_julia()
  fits <- .ac_fits("nonlinear")
  printed <- capture.output(print(fits$on))
  expect_identical(any(grepl("corrected by quadrature", printed)),
    isTRUE(fits$on$laplace$correction$material))
  expect_false(any(grepl("corrected by quadrature", capture.output(print(fits$off)))))
})

test_that("the post-hoc functions do not apply the correction twice", {
  skip_without_julia()
  fits <- .ac_fits("nonlinear")
  on <- fits$on
  expect_error(ctLaplaceCorrect(on), "already corrected")
  expect_error(ctLaplaceCorrect(on, draws = "keep"), "already corrected")
  expect_error(ctLaplaceCorrect(on, draws = "imis", correct_estimate = TRUE),
    "already corrected")

  # The check reports relative to the corrected point: its gap is there, and
  # its correction is the further step, smaller than the one already taken.
  check <- ctLaplaceCheck(on, nodes = 5L)
  expect_identical(check$at, "corrected")
  expect_equal(check$parameters$estimate, as.numeric(on$estimate$raw))
  expect_lt(max(abs(check$parameters$delta)),
    max(abs(on$laplace$correction$delta)))
  expect_true(any(grepl("quadrature-corrected", capture.output(print(check)))))

  # On the uncorrected fit both work as they always did.
  offcheck <- ctLaplaceCheck(fits$off, nodes = 5L)
  expect_identical(offcheck$at, "laplace")
  expect_equal(offcheck$parameters$delta, on$laplace$correction$first_delta,
    tolerance = 1e-6)
  handcorrected <- ctLaplaceCorrect(fits$off, finishsamples = 50)
  expect_s3_class(handcorrected$laplace_correction, "ctLaplaceCorrection")
})

test_that("cross-validation uses the Laplace optimum of a corrected fit", {
  skip_without_julia()
  fits <- .ac_fits("nonlinear")
  expect_equal(.ctLaplaceOptimum(fits$on), fits$on$laplace$correction$laplace_estimate)
  expect_equal(.ctLaplaceOptimum(fits$off), as.numeric(fits$off$estimate$raw))
  loo <- suppressWarnings(suppressMessages(ctLOO(fits$on, folds = 2, cores = 1,
    subjectwise = TRUE)))
  expect_match(loo$scoring, "Laplace optimum")
  # The in-sample side is the Laplace marginal at the Laplace optimum, which is
  # the uncorrected fit's log likelihood, not the corrected fit's.
  expect_equal(loo$insampleLogLik, fits$off$estimate$loglik, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(loo$insampleLogLik, fits$on$estimate$loglik)))
})

test_that("requests that cannot apply are refused by name", {
  # No julia needed.
  expect_error(.ctFitCheckControls(list(laplace_correct = TRUE), "stan"),
    "laplace_correct")
  expect_error(.ctFitCheckControls(list(laplace_correct = FALSE), "stan"), NA)
  expect_error(.ctFitCheckControls(list(laplace_correct = TRUE), "julia"), NA)
  resolve <- function(oc, intoverpop = "laplace", optimize = TRUE,
    intoverstates = TRUE) .ctLaplaceCorrectResolve(oc, intoverpop, optimize,
      intoverstates)
  expect_true(resolve(list()))
  expect_false(resolve(list(laplace_correct = FALSE)))
  expect_error(resolve(list(laplace_correct = "yes")), "TRUE or FALSE")
  expect_error(resolve(list(laplace_correct = NA)), "TRUE or FALSE")
  # Defaulted, it simply does not apply; asked for, it is refused by name.
  expect_false(resolve(list(), intoverpop = "augmented"))
  expect_error(resolve(list(laplace_correct = TRUE), intoverpop = "augmented"),
    "intoverpop='laplace' only")
  expect_false(resolve(list(), optimize = FALSE))
  expect_error(resolve(list(laplace_correct = TRUE), optimize = FALSE),
    "sampled fit")
  expect_error(resolve(list(laplace_correct = TRUE), intoverstates = FALSE),
    "intoverstates=FALSE")
  expect_false(resolve(list(estonly = TRUE)))
  expect_error(resolve(list(laplace_correct = TRUE, estonly = TRUE)), "estonly")
  expect_false(resolve(list(laplace_correct = FALSE, estonly = TRUE)))
})

test_that("refused at ctFit by name, before anything is fitted", {
  skip_without_julia()
  expect_error(ctFit(.ac_linear_data(), .ac_linear_model(), backend = "julia",
    intoverpop = "augmented", optimcontrol = list(laplace_correct = TRUE)),
    "laplace_correct applies to intoverpop='laplace' only")
  estonly <- suppressMessages(ctFit(.ac_linear_data(), .ac_linear_model(),
    backend = "julia", intoverpop = "laplace", optimcontrol = list(estonly = TRUE)))
  expect_identical(estonly$laplace$correction$status, "estonly")
  expect_false(estonly$laplace$correction$applied)
})

test_that("a fit without the record reads as uncorrected", {
  skip_without_julia()
  fits <- .ac_fits("nonlinear")
  # As a fit saved before the default existed: no `$laplace$correction`.
  legacy <- fits$off
  legacy$laplace$correction <- NULL
  expect_false(.ctLaplaceIsCorrected(legacy))
  expect_equal(.ctLaplaceOptimum(legacy), as.numeric(legacy$estimate$raw))
  expect_error(capture.output(print(legacy)), NA)
  check <- ctLaplaceCheck(legacy, nodes = 5L, correction = FALSE)
  expect_identical(check$at, "laplace")
})

test_that("a correction that cannot be evaluated warns and keeps the fit", {
  skip_without_julia()
  fits <- .ac_fits("nonlinear")
  bad <- fits$off
  bad$estimate$raw[1] <- 1e8
  expect_warning(out <- .ctLaplaceAutoCorrect(bad), "Laplace correction skipped")
  expect_identical(out$laplace$correction$status, "failed")
  expect_false(out$laplace$correction$applied)
  expect_equal(out$estimate$raw, bad$estimate$raw)
})

test_that("a gap with no ascent step still reports the quadrature log likelihood", {
  skip_without_julia()
  fits <- .ac_fits("nonlinear")
  # A Hessian with no negative direction: the objective is not concave along
  # any of them, so every direction is dropped and nothing moves -- the shape
  # of AnomAuth's spurious optimum, where the Laplace Hessian was indefinite.
  convex <- fits$off
  convex$uncertainty$hessian <- -convex$uncertainty$hessian
  out <- .ctLaplaceAutoCorrect(convex)
  corr <- out$laplace$correction
  expect_identical(corr$status, "no_gain")
  expect_false(corr$applied)
  expect_equal(corr$dropped_directions, length(convex$estimate$raw))
  expect_equal(out$estimate$raw, convex$estimate$raw)
  expect_identical(out$estimate$loglik_method, "quadrature")
  expect_equal(out$estimate$logposterior,
    .ac_quadrature(out, convex$estimate$raw), tolerance = 1e-8)
  expect_equal(out$estimate$loglik_laplace, convex$estimate$loglik)
  expect_equal(out$estimate$rawposterior, convex$estimate$rawposterior)
  # The print line for this case is keyed on a gap of a nat or more.
  printed <- capture.output(print(out))
  expect_identical(any(grepl("no quadrature step improved", printed)),
    isTRUE(abs(corr$gap_reported) >= 1))
  # The same line when the gap is large, which this fixture's is not.
  big <- out
  big$laplace$correction$gap_reported <- -27.3
  expect_true(any(grepl("exceeds the quadrature value by 27.3",
    capture.output(print(big)))))
})
