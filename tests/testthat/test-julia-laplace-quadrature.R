# The quadrature correction, `optimcontrol$laplace_correct = 'quadrature'` and
# the default (`.ctLaplaceContinue()` in R/ctBackendLaplaceCorrect.R,
# laplace_continuation.jl in the engine), beside the step correction.
#
# What is pinned here: that a model on which Laplace is exact is left alone to
# the bit; that on a nonlinear one the continuation lands at the optimum of the
# EXACT marginal likelihood -- computed here by a dense grid and maximised by
# Nelder-Mead, sharing nothing with the continuation but the model -- closer
# than the step correction does; that what the fit then reports is the
# continuation's own (estimate, quadrature log likelihood, covariance from its
# Hessian, certification); that it is what the default does; that the guard
# reverts; and that each request that cannot apply is refused by name.

# A random effect on an identity-transformed MANIFESTMEANS: Laplace is exact.
.lc_linear_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), MANIFESTMEANS = matrix("mmean"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE
  model
}

.lc_linear_data <- function(nsubjects = 25, nobs = 6) {
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

# One random effect, on a `-log1p_exp` DRIFT; simulated here rather than by
# ctGenerate, whose draw stream moves. The fixture of
# test-julia-laplace-autocorrect.R.
.lc_nonlinear_model <- function() {
  model <- suppressMessages(ctModel(silent = TRUE, type = "ct", CINT = 0,
    MANIFESTMEANS = 0, LAMBDA = matrix(1), T0MEANS = matrix(0),
    DRIFT = "drift|-log1p_exp(-param)|TRUE"))
  model$pars$indvarying <- model$pars$param %in% "drift"
  model
}

.lc_nonlinear_data <- function(seed = 3L, nsubjects = 40L, ntimes = 8L) {
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

.lc_cache <- new.env(parent = emptyenv())
# Every fit of a set shares a seed, so everything before the correction --
# start, optimiser, Hessian, draws -- is the same computation in each.
.lc_fits <- function(which) {
  key <- paste0(which, "_fits")
  if (!exists(key, envir = .lc_cache, inherits = FALSE)) {
    data <- if (identical(which, "linear")) .lc_linear_data() else .lc_nonlinear_data()
    model <- if (identical(which, "linear")) .lc_linear_model() else .lc_nonlinear_model()
    fitwith <- function(correct) {
      set.seed(11)
      suppressMessages(ctFit(data, model, backend = "julia",
        intoverpop = "laplace", cores = 1,
        optimcontrol = list(finishsamples = 100, laplace_correct = correct)))
    }
    fits <- list(off = fitwith(FALSE), quadrature = fitwith("quadrature"))
    if (!identical(which, "linear")) {
      # Asked for by name: the step correction is what the default is
      # compared against.
      fits$step <- fitwith("step")
      fits$default <- fitwith(NULL)
    }
    assign(key, fits, envir = .lc_cache)
  }
  get(key, envir = .lc_cache, inherits = FALSE)
}

# The exact log marginal likelihood of a one-level, one-random-effect Laplace
# objective, by a trapezoid rule over a fixed grid on the standardised effect:
# sum_i log int exp(ll_i(theta + L(theta) u)) phi(u) du, plus the prior. It
# shares nothing with the quadrature but the member shift. 41 points over
# [-7, 7] agree with 161 over [-8, 8] to 1e-10 on this fixture (measured).
.lc_exact <- function(fit, x, npoints = 41L, umax = 7) {
  JuliaConnectoR::juliaEval('
    if !isdefined(Main, :_lc_exact_grid)
      function _lc_exact_grid(laplace, values, umax, npoints)
        C = ContinuousTimeSEM
        theta = collect(Float64, values)
        Ls = C._laplace_popchols(theta, laplace.spec)
        grid = collect(range(-Float64(umax), Float64(umax); length=Int(npoints)))
        h = grid[2] - grid[1]
        logw = [log(h) - u^2 / 2 - log(2pi) / 2 for u in grid]
        logw[1] += log(0.5); logw[end] += log(0.5)
        total = 0.0
        for U in eachindex(laplace.units.members)
          i = laplace.units.members[U][1]
          subject = laplace.objective.subject_objectives[i]
          terms = map(eachindex(grid)) do g
            shifted = C._laplace_member_values(theta, laplace.spec, Ls, [grid[g]],
              laplace.units.offsets[U][1])
            ll = subject(shifted)
            isfinite(ll) ? ll + logw[g] : -Inf
          end
          peak = maximum(terms)
          total += peak + log(sum(t -> exp(t - peak), terms))
        end
        return total + C._ctsem_log_prior(laplace.objective, theta)
      end
    end
    nothing')
  as.numeric(JuliaConnectoR::juliaCall("_lc_exact_grid", .ctJuliaObjective(fit),
    .ctJuliaNumericVector(as.numeric(x)), as.numeric(umax), as.integer(npoints)))
}

test_that("a model on which Laplace is exact is screened and left alone", {
  skip_without_julia()
  fits <- .lc_fits("linear")
  corr <- fits$quadrature$laplace$correction
  expect_identical(corr$method, "quadrature")
  expect_identical(corr$status, "exact")
  expect_false(corr$applied)
  expect_lt(corr$screen, corr$tolerance)
  expect_identical(corr$flagged, 0L)
  # Bit for bit: nothing but the record differs from the uncorrected fit.
  expect_identical(fits$quadrature$estimate$raw, fits$off$estimate$raw)
  expect_identical(fits$quadrature$estimate$loglik, fits$off$estimate$loglik)
  expect_identical(fits$quadrature$estimate$cov, fits$off$estimate$cov)
  expect_identical(fits$quadrature$estimate$rawposterior, fits$off$estimate$rawposterior)
  expect_null(fits$quadrature$estimate$loglik_method)
})

test_that("the default is the quadrature correction, and the fit-time one is the post-hoc one", {
  skip_without_julia()
  fits <- .lc_fits("nonlinear")
  expect_identical(fits$default$laplace$correction$method, "quadrature")
  expect_identical(fits$quadrature$estimate$raw, fits$default$estimate$raw)
  expect_identical(fits$quadrature$estimate$cov, fits$default$estimate$cov)
  expect_identical(fits$quadrature$estimate$rawposterior,
    fits$default$estimate$rawposterior)
  expect_identical(fits$quadrature$estimate$loglik, fits$default$estimate$loglik)
  expect_identical(fits$step$laplace$correction$method, "step")
  # What the fit did is `.ctLaplaceContinue()` applied to the uncorrected fit,
  # exactly, and the same holds for the step.
  posthoc <- .ctLaplaceContinue(fits$off)
  expect_identical(posthoc$estimate$raw, fits$quadrature$estimate$raw)
  expect_identical(posthoc$estimate$loglik, fits$quadrature$estimate$loglik)
  posthoc <- .ctLaplaceAutoCorrect(fits$off)
  expect_identical(posthoc$estimate$raw, fits$step$estimate$raw)
  expect_identical(posthoc$estimate$loglik, fits$step$estimate$loglik)
  # And the post-hoc functions will not apply a correction a second time.
  expect_error(ctLaplaceCorrect(fits$default), "already corrected")
})

test_that("the continuation lands at the exact marginal's optimum, closer than the step", {
  skip_without_julia()
  fits <- .lc_fits("nonlinear")
  cont <- fits$quadrature; step <- fits$step; off <- fits$off
  corr <- cont$laplace$correction
  expect_identical(corr$status, "continued")
  expect_true(corr$applied)
  expect_gte(corr$rounds, 1L)
  expect_identical(corr$continuation, "converged")
  expect_false(corr$guard$fired)
  expect_equal(corr$laplace_estimate, as.numeric(off$estimate$raw))

  # The independent reference: the exact marginal likelihood, maximised
  # without derivatives from the continuation's estimate.
  exact <- function(x) .lc_exact(off, x)
  best <- stats::optim(as.numeric(cont$estimate$raw), function(p) -exact(p),
    method = "Nelder-Mead", control = list(reltol = 1e-11, maxit = 3000))
  best <- stats::optim(best$par, function(p) -exact(p), method = "Nelder-Mead",
    control = list(reltol = 1e-12, maxit = 3000))
  xbest <- best$par
  se <- sqrt(diag(as.matrix(off$estimate$cov)))
  distance <- function(x) max(abs((as.numeric(x) - xbest) / se))
  # Stated tolerance: a tenth of a Laplace standard error in every coordinate.
  # Measured 0.07 for the continuation, 0.23 for the step correction and 1.29
  # for the uncorrected Laplace optimum.
  expect_lt(distance(cont$estimate$raw), 0.1)
  expect_lt(distance(cont$estimate$raw), distance(step$estimate$raw))
  expect_lt(distance(step$estimate$raw), distance(off$estimate$raw))
  # And in exact log likelihood, the order the distances say.
  expect_gt(exact(cont$estimate$raw), exact(step$estimate$raw))
  expect_lt(-best$value - exact(cont$estimate$raw), 0.01)
})

test_that("what a continued fit reports is the continuation's own", {
  skip_without_julia()
  fits <- .lc_fits("nonlinear")
  cont <- fits$quadrature; off <- fits$off
  corr <- cont$laplace$correction
  x <- as.numeric(cont$estimate$raw)
  # The quadrature log likelihood at the estimate, every unit by the rule with
  # its nodes placed there, and the per-subject terms summing to it.
  module <- .ctJuliaModule(cont$model_spec$project)
  quad <- JuliaConnectoR::juliaGet(module$ctsem_laplace_quadrature(
    .ctJuliaObjective(cont), .ctJuliaNumericVector(x), nodes = 5L))$value
  expect_identical(cont$estimate$loglik_method, "quadrature")
  expect_equal(cont$estimate$logposterior, quad, tolerance = 1e-8)
  expect_equal(cont$estimate$loglik, corr$loglik_quadrature)
  expect_equal(sum(cont$estimate$subject_loglik), cont$estimate$loglik, tolerance = 1e-8)
  expect_equal(cont$estimate$loglik_laplace, off$estimate$loglik)
  # The covariance and draws come from the continuation's Hessian, evaluated at
  # its estimate, not from the Laplace curvature at the Laplace optimum.
  expect_equal(as.numeric(cont$uncertainty$evaluated_at), x)
  expect_equal(cont$uncertainty$hessian, corr$hessian)
  expect_false(isTRUE(all.equal(as.matrix(cont$estimate$cov),
    as.matrix(off$estimate$cov))))
  # Values only: the covariance carries its construction's diagnostics as an
  # attribute.
  expect_equal(as.numeric(cont$estimate$cov),
    as.numeric(solve(-(corr$hessian + t(corr$hessian)) / 2)), tolerance = 1e-6)
  expect_lt(max(abs(colMeans(cont$estimate$rawposterior) - x) /
    as.numeric(cont$estimate$se)), 0.5)
  expect_identical(corr$draws, "redrawn")
  # Certified on the objective it maximises.
  expect_true(isTRUE(corr$certification$certified))
  expect_identical(cont$uncertainty$certification$status, corr$certification$status)
  expect_true(isTRUE(cont$optim$converged))
  # The accounting a reader needs.
  expect_true(all(c("values", "gradients", "placements") %in% names(corr$evaluations)))
  expect_gt(corr$flagged, 0L)
  expect_lte(corr$residual, corr$tolerance)
  expect_true(is.data.frame(corr$trace))
  # The Laplace optimum stays what cross-validation reads.
  expect_equal(.ctLaplaceOptimum(cont), as.numeric(off$estimate$raw))
  expect_true(.ctLaplaceIsCorrected(cont))
  # And print says what ran.
  printed <- capture.output(print(cont))
  expect_identical(any(grepl("continued on the quadrature objective", printed)),
    isTRUE(corr$material))
})

test_that("ctLaplaceCheck agrees with a continued fit and refines by continuing", {
  skip_without_julia()
  fits <- .lc_fits("nonlinear")
  cont <- fits$quadrature; off <- fits$off
  # At the continuation's estimate the further first-order step is small but
  # not zero, and it should not be: the check steps towards the maximum of the
  # quadrature VALUE with its nodes re-placed at every point, and the
  # continuation stops at the zero of the quadrature estimate of the score,
  # which on this fixture is the nearer of the two to the exact optimum (0.07
  # against 0.14 standard errors). Measured 0.057.
  check <- ctLaplaceCheck(cont, nodes = 5L)
  expect_identical(check$at, "corrected")
  expect_lt(max(abs(check$parameters$delta_se), na.rm = TRUE), 0.1)
  # The two reports count the same directions.
  expect_identical(check$dropped_directions, as.integer(cont$identifiability$nweak))
  # refine = TRUE on the uncorrected fit is the same continuation, from the
  # same start and with the same curvature, so it lands where the fit did.
  refined <- ctLaplaceCheck(off, nodes = 5L, refine = TRUE)
  expect_identical(refined$refined_status, "converged")
  expect_true(refined$refined_converged)
  expect_equal(refined$refined, as.numeric(cont$estimate$raw), tolerance = 1e-8)
})

test_that("the guard reverts a continuation that moves the objective too far", {
  skip_without_julia()
  fits <- .lc_fits("nonlinear")
  off <- fits$off
  control <- .ctLaplaceContinueDefaults
  control$guard <- 1e-6
  control$guard_per_subject <- 0
  expect_warning(out <- .ctLaplaceContinue(off, control = control),
    "keeps the Laplace optimum")
  corr <- out$laplace$correction
  expect_identical(corr$status, "reverted")
  expect_false(corr$applied)
  expect_true(corr$guard$fired)
  expect_identical(out$estimate$raw, off$estimate$raw)
  expect_identical(out$estimate$cov, off$estimate$cov)
  expect_false(isTRUE(all.equal(corr$rejected_estimate, as.numeric(off$estimate$raw))))
  # The quadrature log likelihood at the Laplace optimum is what is reported.
  expect_identical(out$estimate$loglik_method, "quadrature")
})

test_that("requests that cannot apply are refused by name", {
  # No julia needed.
  resolve <- function(oc, intoverpop = "laplace", optimize = TRUE,
    intoverstates = TRUE) .ctLaplaceCorrectResolve(oc, intoverpop, optimize,
      intoverstates)
  expect_identical(resolve(list()), "quadrature")
  expect_identical(resolve(list(laplace_correct = TRUE)), "quadrature")
  expect_identical(resolve(list(laplace_correct = "step")), "step")
  expect_identical(resolve(list(laplace_correct = "quadrature")), "quadrature")
  expect_false(resolve(list(laplace_correct = FALSE)))
  expect_error(resolve(list(laplace_correct = "yes")), "'step' or 'quadrature'")
  expect_error(resolve(list(laplace_correct = "continue")), "'step' or 'quadrature'")
  expect_error(resolve(list(laplace_correct = c("step", "quadrature"))),
    "'step' or 'quadrature'")
  expect_error(resolve(list(laplace_correct = "quadrature"), intoverpop = "augmented"),
    "intoverpop='laplace' only")
  expect_error(resolve(list(laplace_correct = "quadrature"), optimize = FALSE),
    "sampled fit")
  expect_error(resolve(list(laplace_correct = "quadrature", estonly = TRUE)), "estonly")
  expect_error(.ctFitCheckControls(list(laplace_correct = "quadrature"), "stan"),
    "laplace_correct")
})

# The gated-gaps A14 config (review/LAPLACE-gated-gaps-2026-09-24.md): 40
# subjects of 6 waves, a random -log1p_exp DRIFT beside random T0MEANS and
# CINT, so three effects per subject and the soft-direction rule. Its data
# generator, copied from the gaps job's (dev/optimbench/cells.R, gg_genA).
.lc_a14_data <- function(seed = 2L, nsub = 40L, ntimes = 6L, mu = 3, sdp = 2.5) {
  set.seed(seed)
  baseline <- stats::rnorm(nsub, 2, 2)
  start <- stats::rnorm(nsub, baseline / 2, 1)
  raw <- stats::rnorm(nsub, mu + (baseline - 2) / 2, sdp)
  drift <- -log1p(exp(-raw))
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    a <- drift[i]; decay <- exp(a)
    intercept <- (baseline[i] / a) * (decay - 1)
    innovation <- sqrt(0.25 * (exp(2 * a) - 1) / (2 * a))
    latent <- numeric(ntimes); latent[1] <- start[i]
    for (t in seq_len(ntimes - 1L)) latent[t + 1L] <- decay * latent[t] +
      intercept + stats::rnorm(1, 0, innovation)
    data.frame(id = i, time = seq_len(ntimes) - 1L,
      Y1 = latent + stats::rnorm(ntimes, 0, 0.5))
  }))
}

test_that("a continuation whose fixed-node model misleads it does not walk downhill", {
  skip_without_julia()
  # Slow because the failure needs the soft rule, so three effects a subject,
  # and this is the data it was measured on. The fixture above has two effects,
  # where the product rule applies, and passed under the old acceptance. A
  # smaller version of this config was not tried.
  skip_unless_slow("the A14 continuation (a 40-subject, three-effect Laplace fit)")
  model <- suppressMessages(ctModel(silent = TRUE, type = "ct", CINT = "cint",
    MANIFESTMEANS = 0, LAMBDA = matrix(1), DRIFT = "drift|-log1p_exp(-param)|TRUE"))
  set.seed(1)
  fit <- suppressWarnings(suppressMessages(ctFit(.lc_a14_data(), model,
    backend = "julia", intoverpop = "laplace", cores = 1,
    optimcontrol = list(finishsamples = 100))))
  corr <- fit$laplace$correction
  tol <- .ctLaplaceContinueDefaults$value_tol
  expect_true(corr$status %in% c("continued", "no_gain"))
  # Kept rounds used to be those that lowered the fixed-point residual alone.
  # Here that gave up 0.9 nats of re-placed value on the way to the fixed
  # point, which is 0.84 exact nats below where the rounds now stop (-398.81
  # against -397.97; Laplace -399.23); with the one-node stiff complement the
  # soft rule had before, the estimate ended 4.9 exact nats below Laplace. Now
  # a kept round may not lower the value, and a run that ends lower than it
  # began keeps the Laplace optimum.
  expect_true(all(corr$trace$gain[corr$trace$kept] >= -tol))
  if (identical(corr$status, "continued")) {
    expect_gte(corr$guard$change, -tol)
  } else {
    expect_equal(fit$estimate$raw, corr$laplace_estimate)
  }
  # The rounds stop short of the fixed point here. The continuation's Hessian
  # at such a point had two directions of positive curvature, which left ten
  # parameters with no interval and the fit `converged: FALSE`; so a
  # continuation that stopped short keeps the fit's own uncertainty, at the
  # Laplace optimum, and recentres its draws.
  if (isTRUE(corr$applied) && !identical(corr$continuation, "converged")) {
    expect_identical(corr$hessian_at, "laplace_estimate")
    expect_identical(corr$draws, "recentred")
    expect_equal(as.numeric(fit$uncertainty$evaluated_at),
      as.numeric(corr$laplace_estimate))
  }
})
