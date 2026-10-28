# The quadrature correction, `optimcontrol$laplace_correct = 'quadrature'`,
# and the default's check (`.ctLaplaceContinue()` in
# R/ctBackendLaplaceCorrect.R, laplace_continuation.jl in the engine), beside
# the step correction.
#
# What is pinned here: that a model on which Laplace is exact is left alone to
# the bit; that on a nonlinear one the continuation lands at the optimum of the
# EXACT marginal likelihood -- computed here by a dense grid and maximised by
# Nelder-Mead, sharing nothing with the continuation but the model -- closer
# than the step correction does; that what the fit then reports is the
# continuation's own (estimate, quadrature log likelihood, covariance from its
# Hessian, certification); that the default checks, warns and leaves the
# estimate, and that the correction its warning names is the fit-time one;
# that the guard reverts; and that each request that cannot apply is refused
# by name.

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
    # 'checked' by name: what these tests read is the draws a corrected fit
    # takes on the hybrid, which the default's curvature does not take.
    # The warnings are kept on the fit, for the tests that read them: the
    # default's check warns where it finds a move.
    fitwith <- function(correct, uncertainty = "checked") {
      set.seed(11)
      warned <- character()
      fit <- withCallingHandlers(suppressMessages(ctFit(data, model,
        backend = "julia", intoverpop = "laplace", cores = 1,
        optimcontrol = list(finishsamples = 100, laplace_correct = correct,
          uncertainty = uncertainty))),
        warning = function(w) {
          warned <<- c(warned, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      attr(fit, "warned") <- warned
      fit
    }
    fits <- list(off = fitwith(FALSE), quadrature = fitwith("quadrature"))
    if (!identical(which, "linear")) {
      # Asked for by name: the step correction is what the quadrature
      # correction is compared against.
      fits$step <- fitwith("step")
      fits$default <- fitwith(NULL)
      # Every default: the curvature's normal, and the shift check alone.
      fits$development <- fitwith(NULL, uncertainty = NULL)
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

test_that("the default checks and warns, and the correction it names is the fit-time one", {
  skip_without_julia()
  fits <- .lc_fits("nonlinear")
  default <- fits$default
  corr <- default$laplace$correction
  # Correcting Laplace moves the population sd here, resolved past the
  # check's noise, so the default says so, names the corrections, and leaves
  # the estimate, its log likelihood and its draws as the uncorrected fit's.
  expect_identical(corr$method, "quadrature")
  expect_identical(corr$status, "flagged")
  expect_false(corr$applied)
  expect_identical(corr$check$status, "flagged")
  expect_identical(corr$check$wide, 0L)
  expect_false(grepl("hours", corr$message, fixed = TRUE))
  expect_gt(corr$check$carrying, 0L)
  expect_length(corr$check$carriers, corr$check$carrying)
  expect_identical(default$estimate$raw, fits$off$estimate$raw)
  expect_identical(default$estimate$loglik, fits$off$estimate$loglik)
  expect_identical(default$estimate$rawposterior, fits$off$estimate$rawposterior)
  expect_null(default$estimate$loglik_method)
  warned <- attr(default, "warned")
  expect_true(any(grepl("likelihood = 'quadrature'", warned, fixed = TRUE)))
  expect_true(any(grepl("intoverpop = 'none'", warned, fixed = TRUE)))
  expect_true(any(grepl(corr$message, capture.output(print(default)), fixed = TRUE)))
  # The correction it names is the fit-time one: the rounds on the units its
  # check named, from the same estimate, then the hybrid's own curvature.
  corrected <- suppressMessages(ctFitUncertainty(default, "hessian",
    control = list(likelihood = "quadrature")))
  expect_identical(corrected$laplace$correction$status, "continued")
  expect_equal(corrected$estimate$raw, fits$quadrature$estimate$raw, tolerance = 1e-6)
  expect_equal(corrected$estimate$loglik, fits$quadrature$estimate$loglik,
    tolerance = 1e-6)
  expect_identical(corrected$uncertainty$settings$method, "hessian")
  expect_identical(corrected$uncertainty$settings$draws, "normal")
  expect_identical(fits$step$laplace$correction$method, "step")
  # The Laplace run converged (no unit's term jumped, `fit$laplace$jump`), so
  # what the fit did is `.ctLaplaceContinue()` applied to the uncorrected fit,
  # exactly, and the same holds for the step.
  expect_null(fits$quadrature$laplace$jump)
  posthoc <- .ctLaplaceContinue(fits$off)
  expect_identical(posthoc$estimate$raw, fits$quadrature$estimate$raw)
  expect_identical(posthoc$estimate$loglik, fits$quadrature$estimate$loglik)
  posthoc <- .ctLaplaceAutoCorrect(fits$off)
  expect_identical(posthoc$estimate$raw, fits$step$estimate$raw)
  expect_identical(posthoc$estimate$loglik, fits$step$estimate$loglik)
  # And sampling against the quadrature will not apply it a second time.
  expect_error(ctFitUncertainty(fits$quadrature, "is", control = list(
    likelihood = "quadrature", correctEstimate = TRUE)), "already corrected")
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
  # Converged means the gain the fixed-node model still predicts is under the
  # bar the rounds stopped on, in nats: the certification tolerance, while
  # `stop_gain` is off.
  expect_lt(corr$stationarity, corr$stop_gain)
  expect_identical(corr$stop_gain, .ctBackendGapTolerance(off))
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
  # The Laplace log likelihood kept is the Laplace objective's at the point the
  # correction started from, the climb's end.
  lap <- JuliaConnectoR::juliaGet(module$ctsem_laplace_evaluate(
    .ctJuliaObjective(cont), .ctJuliaNumericVector(corr$laplace_estimate),
    gradient = FALSE))
  expect_equal(cont$estimate$loglik_laplace, sum(lap$unit_loglik), tolerance = 1e-8)
  # The covariance and draws come from the continuation's Hessian, evaluated at
  # its estimate -- or where the finish's steps started, within the reuse
  # distance of it -- not from the Laplace curvature at the Laplace optimum.
  expect_lte(.ctBackendHessianDistance(corr$hessian,
    cont$uncertainty$evaluated_at, x), .ctBackendHessianReuse())
  expect_identical(corr$newton$steps > 0L,
    !isTRUE(all.equal(as.numeric(cont$uncertainty$evaluated_at), x, tolerance = 0)))
  expect_equal(cont$uncertainty$hessian, corr$hessian)
  expect_false(isTRUE(all.equal(as.matrix(cont$estimate$cov),
    as.matrix(off$estimate$cov))))
  # The normal is the continuation's Hessian's inverse. 'checked', asked for
  # here, reports it unless its check finds it wrong along some direction; then it
  # importance-samples the continuation's quadrature objective from it and
  # reports the weighted draws' covariance, the normal's construction kept as
  # the proposal's (`details$proposal_covariance`). On this fixture the check
  # finds the normal wrong (a slope of -4.4 along one direction), so the
  # covariance is the sampler's, and the proposal is tied to the Hessian by
  # the information's smallest eigenvalue, which its record keeps. Values only:
  # the covariance carries its construction's diagnostics as an attribute.
  information <- -(corr$hessian + t(corr$hessian)) / 2
  details <- cont$uncertainty$details
  # The method asked for, drawn once, on the continuation's hybrid: the pass
  # before the continuation took only the curvature's normal.
  expect_identical(cont$uncertainty$settings$method, "checked")
  if (isTRUE(details$importance_sampling$used)) {
    expect_true(isTRUE(details$importance_sampling$weighted))
    # The hybrid's density carries the gradient batches the Stein check reads.
    expect_false(is.null(details$importance_sampling$stein))
    expect_identical(details$proposal_covariance$method, "solve")
    expect_false(isTRUE(details$proposal_covariance$usedNullProjection))
    expect_equal(details$proposal_covariance$minInfoEigenOriginal,
      min(eigen(information, symmetric = TRUE, only.values = TRUE)$values),
      tolerance = 1e-6)
  } else {
    expect_equal(as.numeric(cont$estimate$cov), as.numeric(solve(information)),
      tolerance = 1e-6)
  }
  # The draws are the continuation's posterior's, not the Laplace optimum's:
  # in every coordinate where the two estimates separate by more than the
  # draws' error, the draws' mean is the continued estimate's side. Not within
  # a fixed share of a standard error of the estimate: 'checked' centres them
  # on the posterior mean (b063aaaf), which on this fixture sits 0.5-0.7 se
  # from the mode in two coordinates where the two estimates agree.
  se <- as.numeric(cont$estimate$se)
  centre <- colMeans(cont$estimate$rawposterior)
  laplace <- as.numeric(off$estimate$raw)
  apart <- abs(x - laplace) / se > 0.3
  expect_true(any(apart))
  expect_true(all(abs(centre - x)[apart] < abs(centre - laplace)[apart]))
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
  # What cross-validation reads is the point the correction started from.
  expect_equal(.ctLaplaceOptimum(cont), as.numeric(corr$laplace_estimate))
  expect_true(.ctLaplaceIsCorrected(cont))
  # And print says what ran.
  printed <- capture.output(print(cont))
  expect_identical(any(grepl("continued on the quadrature objective", printed)),
    isTRUE(corr$material))
})

test_that("a default fit keeps the Laplace curvature's normal, beside both checks", {
  skip_without_julia()
  fits <- .lc_fits("nonlinear")
  dev <- fits$development
  corr <- dev$laplace$correction
  # Checked as the fit asking for 'checked' was, on the same curvature, and
  # left where it is.
  expect_identical(corr$status, "flagged")
  expect_equal(corr$check$shift, fits$default$laplace$correction$check$shift)
  expect_equal(dev$estimate$raw, fits$off$estimate$raw)
  expect_identical(dev$uncertainty$settings$method, "hessian")
  # The normal is the Laplace curvature's at the estimate, and the normal
  # check, made beside it, is kept and changes nothing.
  hessian <- dev$uncertainty$hessian
  information <- -(hessian + t(hessian)) / 2
  expect_equal(as.numeric(dev$estimate$cov), as.numeric(solve(information)),
    tolerance = 1e-6)
  check <- dev$uncertainty$details$normalcheck
  expect_false(is.null(check))
  expect_null(dev$uncertainty$details$importance_sampling)
  # And summary() says what it found, where it found something.
  named <- if (is.matrix(check$directions)) ncol(check$directions) else 0L
  note <- suppressMessages(summary(dev, parmatrices = FALSE))$uncertaintyNote
  expect_identical(grepl("ctFitUncertainty(fit, 'sample')", note, fixed = TRUE),
    named > 0L)
})

test_that("the hybrid's density places the flagged units' rules at each point", {
  skip_without_julia()
  fit <- .lc_fits("nonlinear")$quadrature
  module <- .ctJuliaModule(fit$model_spec$project)
  x <- as.numeric(fit$estimate$raw)
  cont <- do.call(module$ctsem_laplace_continuation,
    c(list(.ctJuliaObjective(fit, hybrid = TRUE), .ctJuliaNumericVector(x)),
      .ctLaplaceContinuationArgs(.ctLaplaceContinueDefaults)))
  info <- .ctJuliaGet(module$ctsem_laplace_continuation_info(cont))
  flagged <- as.integer(info$flagged)[seq_len(info$nflagged)]
  expect_gt(length(flagged), 0L)
  frozen <- function(z) .ctJuliaGet(module$ctsem_laplace_continuation_evaluate(cont,
    .ctJuliaNumericVector(z), gradient = TRUE))$value
  before <- frozen(x)
  lpg <- .ctLaplaceContinueLpg(module, cont)(fit)
  away <- x + 2 * as.numeric(fit$estimate$se) * rep_len(c(1, -1), length(x))
  value <- lpg(away)
  # The same as a continuation placed there from nothing, with the same units
  # flagged and no screen: the hybrid's value less its frozen draws' error.
  there <- .ctJuliaGet(module$ctsem_laplace_continuation_info(do.call(
    module$ctsem_laplace_continuation, c(list(.ctJuliaObjective(fit, hybrid = TRUE),
      .ctJuliaNumericVector(away)), .ctLaplaceContinuationArgs(
        .ctLaplaceContinueDefaults, force = flagged, screen = FALSE)))))
  expect_equal(as.numeric(value), as.numeric(there$value) - as.numeric(there$error_draws),
    tolerance = 1e-6)
  expect_length(attr(value, "gradient"), length(x))
  expect_true(all(is.finite(attr(value, "gradient"))))
  # A batch is the points one at a time, values with or without gradients.
  both <- rbind(away, x)
  batch <- attr(lpg, "batch")(both)
  expect_equal(batch[1L], as.numeric(value), tolerance = 1e-8)
  expect_equal(attr(lpg, "gradbatch")(both)$value, batch, tolerance = 1e-8)
  # And the continuation is left placed where it was.
  after <- .ctJuliaGet(module$ctsem_laplace_continuation_info(cont))
  expect_equal(as.numeric(after$centre), as.numeric(info$centre))
  expect_identical(after$recentres, info$recentres)
  expect_equal(frozen(x), before, tolerance = 1e-10)
})

test_that("drawing a corrected fit on its hybrid leaves it corrected", {
  skip_without_julia()
  skip_if_not_installed('mvtnorm'); skip_if_not_installed('diagis')
  fit <- .lc_fits("nonlinear")$quadrature
  set.seed(2)
  redrawn <- suppressWarnings(ctFitUncertainty(fit, "is", finishsamples = 20,
    control = list(likelihood = "quadrature", isitersize = 20, imisMaxIter = 1L,
      isESS = 10)))
  # Drawn on the hybrid, placed where the fit stands: before, the run that
  # moved nothing returned before drawing, and its record replaced the fit's.
  expect_false(is.null(redrawn$uncertainty$details$laplace_correction$error))
  expect_identical(redrawn$laplace$correction, fit$laplace$correction)
  expect_equal(as.numeric(redrawn$estimate$raw), as.numeric(fit$estimate$raw))
  expect_identical(redrawn$uncertainty$settings$method, "is")
})

test_that("the continuation closes on its own Hessian what its rounds' rule left", {
  skip_without_julia()
  # The rounds stop on a residual in the Laplace fit's metric, and the
  # certification measures the gap on the continuation's Hessian; where that
  # Hessian says more is left, `.ctLaplaceContinueFinish()` takes the steps it
  # predicts. Here from points a little off the continued estimate, with the
  # estimate's Hessian standing in for theirs, as it does within the reuse
  # distance.
  fit <- .lc_fits("nonlinear")$quadrature
  corr <- fit$laplace$correction
  module <- .ctJuliaModule(fit$model_spec$project)
  control <- .ctLaplaceContinueDefaults
  x <- as.numeric(fit$estimate$raw)
  hc <- corr$hessian
  se <- sqrt(diag(solve(-(hc + t(hc)) / 2)))
  signs <- rep_len(c(1, -1), length(x))
  placed_at <- function(z) {
    cont <- module$ctsem_laplace_continuation(.ctJuliaObjective(fit, hybrid = TRUE),
      .ctJuliaNumericVector(z), nodes = as.integer(control$nodes),
      tolerance = as.numeric(control$tolerance),
      product_maxdim = as.integer(control$product_maxdim),
      soft_tau = as.numeric(control$soft_tau),
      soft_maxdirs = as.integer(control$soft_maxdirs),
      maxdim = as.integer(control$maxdim))
    list(cont = cont, final = .ctJuliaGet(module$ctsem_laplace_continuation_evaluate(
      cont, .ctJuliaNumericVector(z), gradient = TRUE)),
      value = as.numeric(.ctJuliaGet(
        module$ctsem_laplace_continuation_info(cont))$quadrature))
  }

  near <- x + 0.004 * signs * se
  at <- placed_at(near)
  out <- .ctLaplaceContinueFinish(module, at$cont, near, hc, at$final,
    value = at$value, bar = 1e-12, value_tol = control$value_tol)
  expect_gte(out$steps, 1L)
  expect_gt(out$gap_before, 1e-6)
  expect_lt(out$gap_after, out$gap_before / 100)
  expect_gte(out$value, at$value)
  # Back where the rounds converged, to the precision their bar left it at.
  expect_lt(max(abs(out$x - x) / se), 0.003)
  # What it returns is the continuation placed at its point.
  expect_equal(as.numeric(.ctJuliaGet(
    module$ctsem_laplace_continuation_info(at$cont))$quadrature), out$value)

  # A step farther than the Hessian describes is not taken.
  far <- x + 0.2 * signs * se
  at <- placed_at(far)
  out <- .ctLaplaceContinueFinish(module, at$cont, far, hc, at$final,
    value = at$value, bar = 1e-12, value_tol = control$value_tol)
  expect_identical(out$steps, 0L)
  expect_identical(out$x, far)
  expect_null(out$info)
})

test_that("the quadrature correction reports progress at the default verbosity", {
  skip_without_julia()
  fits <- .lc_fits("nonlinear")
  # `.ctProgressConsole()` mocked TRUE stands in for a real console this test
  # cannot open: reporting is meant to appear whenever something is watching,
  # by default, with nothing asked for on `verbose` -- unmocked, a testthat
  # run is non-interactive and would report nothing, which is correct there
  # and would prove nothing here.
  testthat::local_mocked_bindings(.ctProgressConsole = function() TRUE)
  seen <- character()
  withCallingHandlers(
    invisible(.ctLaplaceContinue(fits$off)),
    message = function(m) {
      seen <<- c(seen, conditionMessage(m)); invokeRestart("muffleMessage")
    })
  # Screen, at least one round (this fixture converges in >= 1, per the
  # continuation test above) and the final Hessian -- the three phases that
  # went a combined 22 minutes without a word on Charles's ordinal fit
  # (2026-09-28) because each was gated on `verbose > 0` alone.
  expect_true(any(grepl("Laplace continuation screen", seen, fixed = TRUE)))
  expect_true(any(grepl("Laplace continuation round", seen, fixed = TRUE)))
  expect_true(any(grepl("Laplace continuation hessian", seen, fixed = TRUE)))
})

test_that("a certification that forms its own Hessian reports progress through the shared sink", {
  skip_without_julia()
  # Fit-free: a certification that forms its own Hessian is `ctsem_endgame`
  # (R/ctBackendOptimGap.R's `.ctBackendEndgameAt()`), reached here directly
  # with no optimisation run at all -- only an objective, built from an
  # already-fitted spec's Laplace layout, and a point to evaluate at.
  fits <- .lc_fits("nonlinear")
  spec <- structure(fits$off$model_spec, class = c("ctJuliaModel", "ctFitModel"))
  module <- ctsem:::.ctJuliaModule(spec$project)
  objective <- ctsem:::.ctJuliaObjective(spec)
  est <- as.numeric(fits$off$estimate$raw)
  seen <- character()
  capture <- function(text, kind = "update") seen <<- c(seen, text)
  # `progress_every` set far below any real cadence: this fixture's Hessian is
  # only a handful of columns, fast enough on an already-compiled shape that
  # the production cadence (0.4s/5s) might never be crossed inside a test's
  # patience. That is not what is under test here -- the finish already
  # shares this reporter with the main optimiser's own line, which uses the
  # production cadence -- what is under test is that the callback reaches the
  # sink at all, with the right counts.
  out <- JuliaConnectoR::juliaGet(module$ctsem_endgame(objective,
    ctsem:::.ctJuliaNumericVector(est), progress = TRUE, progress_sink = capture,
    progress_overwrite = TRUE, progress_label = "certify",
    progress_every = 1e-9))
  npar <- length(est)
  expect_true(is.matrix(out$hessian))
  expect_equal(dim(out$hessian), c(npar, npar))
  expect_true(any(grepl("hessian", seen, fixed = TRUE)))
  expect_true(any(grepl(paste0("of ", 2L * npar, " gradients"), seen)))
  # The closing line, once something was actually shown.
  expect_true(any(grepl("hessian formed", seen, fixed = TRUE)))
})

test_that("ctsem_hessian_progress reports the same way, for the caller with no finish to borrow a reporter from", {
  skip_without_julia()
  # `.ctBackendHessian()` (R/ctBackendUncertainty.R) is the standalone
  # `ctFitUncertainty()`/certification Hessian call: reached with no running
  # optimiser and no `ctsem_endgame` certification either. Before this it hand
  # -built an R-side rate limiter and a copy of the "hessian K of N gradients"
  # text `ctsem_endgame` already produces; now it asks for the engine's own
  # reporter the same way `ctsem_endgame` does, through `ctsem_hessian_progress`.
  fits <- .lc_fits("nonlinear")
  spec <- structure(fits$off$model_spec, class = c("ctJuliaModel", "ctFitModel"))
  module <- ctsem:::.ctJuliaModule(spec$project)
  objective <- ctsem:::.ctJuliaObjective(spec)
  est <- as.numeric(fits$off$estimate$raw)
  npar <- length(est)

  seen <- character()
  capture <- function(text, kind = "update") seen <<- c(seen, text)
  # `ctsem_hessian`/`ctsem_hessian_forward` return a bare `Matrix{Float64}`,
  # which JuliaConnectoR marshals straight to an R matrix rather than a proxy
  # -- `.ctBackendJuliaValue()` is what the rest of the package calls a value
  # that may or may not need `juliaGet()` through, and is what
  # `.ctBackendHessian()` itself uses below.
  result <- ctsem:::.ctBackendJuliaValue(module$ctsem_hessian_progress(objective,
    ctsem:::.ctJuliaNumericVector(est), progress = TRUE, progress_sink = capture,
    progress_overwrite = TRUE, progress_label = "hessian", progress_every = 1e-9))
  hessian <- matrix(as.numeric(result), npar, npar)
  expect_true(all(is.finite(hessian)))
  expect_true(any(grepl(paste0("of ", 2L * npar, " gradients"), seen)))

  # Reporting is a side channel, not a second computation: the same call with
  # no reporter at all gives the same matrix. `expect_equal()`, not
  # `expect_identical()`: this call is not pinned to one chunk, and two
  # separate invocations summing per-unit contributions over however many
  # threads the session has can differ in the last bit or two from
  # floating-point addition's own non-associativity -- observed here as
  # "equal but not identical", on this direct call, which sets no thread or
  # cores ceiling of its own. A real second computation would move far more
  # than that.
  bare <- matrix(as.numeric(ctsem:::.ctBackendJuliaValue(
    module$ctsem_hessian_progress(objective, ctsem:::.ctJuliaNumericVector(est)))),
    npar, npar)
  expect_equal(hessian, bare, tolerance = 1e-10)

  # And `.ctBackendHessian()` reaches this same function and the same numbers.
  # `estonly = TRUE` leaves the fit with no stored Hessian at all
  # (`fit$uncertainty` is NULL), so `.ctBackendStoredHessian()` cannot short
  # -circuit the call and this actually exercises the computation, not a
  # cache hit.
  estonly <- suppressMessages(ctFit(.lc_nonlinear_data(), .lc_nonlinear_model(),
    backend = "julia", intoverpop = "laplace", cores = 1,
    optimcontrol = list(estonly = TRUE)))
  expect_null(estonly$uncertainty)
  direct <- suppressMessages(ctsem:::.ctBackendHessian(estonly,
    as.numeric(estonly$estimate$raw)))
  expect_true(is.matrix(direct))
  expect_equal(dim(direct), c(npar, npar))
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

test_that("a correction predicted to gain less than skip_gain is skipped, and says so", {
  skip_without_julia()
  fits <- .lc_fits("nonlinear")
  # This fixture's correction predicts about 0.2 nats from its Laplace
  # optimum, so at the default bar it runs; a bar of a nat skips it.
  expect_gte(fits$quadrature$laplace$correction$predicted_gain,
    .ctLaplaceContinueDefaults$skip_gain)
  high <- utils::modifyList(.ctLaplaceContinueDefaults, list(skip_gain = 1))
  out <- .ctLaplaceContinue(fits$off, control = high)
  corr <- out$laplace$correction
  expect_identical(corr$status, "skipped")
  expect_false(corr$applied)
  expect_equal(corr$predicted_gain,
    fits$quadrature$laplace$correction$predicted_gain)
  expect_null(corr$trace)
  # Nothing moved: the estimate and its uncertainty are the Laplace fit's, and
  # the log likelihood reported is the quadrature one there.
  expect_identical(out$estimate$raw, fits$off$estimate$raw)
  expect_identical(out$estimate$cov, fits$off$estimate$cov)
  expect_identical(out$estimate$loglik_method, "quadrature")
  expect_equal(out$estimate$loglik, corr$loglik_quadrature)
  expect_true(any(grepl("the quadrature correction predicted a gain",
    capture.output(print(out)))))
})

test_that("a unit wider than the cap keeps its Laplace term, and the fit says so", {
  skip_without_julia()
  # One effect a subject on a nonlinear drift; a cap of none leaves every
  # subject unscored by the screen, with draws off or kept, by default, to the
  # units a rescue or the shift check named. With them for every unit, a unit
  # of one block past the cap is scored by frozen draws instead.
  fits <- .lc_fits("nonlinear")
  off <- fits$off
  narrow <- utils::modifyList(.ctLaplaceContinueDefaults, list(maxdim = 0L,
    draws = 0L, check = FALSE))
  expect_message(out <- .ctLaplaceContinue(off, control = narrow),
    "kept the Laplace term, having more than 0 random effects")
  corr <- out$laplace$correction
  expect_identical(corr$status, "too_wide")
  expect_false(corr$applied)
  expect_identical(corr$wide, corr$units)
  expect_identical(out$estimate$raw, off$estimate$raw)
  expect_identical(out$estimate$loglik, off$estimate$loglik)
  expect_true(any(grepl("not corrected: every subject or group has more than 0",
    capture.output(print(out)))))
  # Under the default cap nothing here is wide.
  expect_identical(fits$quadrature$laplace$correction$wide, 0L)
  # Nor past it, with draws: each subject is one block, scored by them. Read
  # after the screen, before any round.
  drawn <- suppressMessages(.ctLaplaceContinue(off, control = utils::modifyList(
    narrow, list(draws = 256L, draws_scope = "all", skip_gain = Inf))))
  expect_identical(drawn$laplace$correction$wide, 0L)
  expect_true(drawn$laplace$correction$status %in% c("skipped", "exact"))
  # The shift check takes no notice of the cap: it scores every unit on the
  # same frozen draws, and reads the same step as the default fit's did.
  checked <- suppressWarnings(suppressMessages(.ctLaplaceContinue(off,
    control = utils::modifyList(.ctLaplaceContinueDefaults,
      list(maxdim = 0L, skip_gain = Inf)))))
  expect_equal(checked$laplace$correction$check$shift,
    fits$quadrature$laplace$correction$check$shift)
  # On a Gaussian integrand Laplace is exact at any width: no check, no
  # screen, and nothing reported as left on Laplace -- with the effect a
  # coordinate as much as with the subject filter integrating it.
  linear <- .lc_fits("linear")$off
  layer <- linear
  layer$model_spec$laplace$hybrid <- NULL
  for (f in list(linear, layer)) {
    expect_no_message(exact <- .ctLaplaceContinue(f, control = narrow))
    expect_identical(exact$laplace$correction$status, "exact")
    expect_identical(exact$laplace$correction$wide, 0L)
    expect_null(exact$laplace$correction$check)
    expect_identical(exact$estimate$raw, f$estimate$raw)
  }
})

test_that("requests that cannot apply are refused by name", {
  # No julia needed.
  resolve <- function(oc, intoverpop = "laplace", optimize = TRUE,
    intoverstates = TRUE) .ctLaplaceCorrectResolve(oc, intoverpop, optimize,
      intoverstates)
  expect_identical(resolve(list()), "check")
  expect_identical(resolve(list(laplace_correct = "check")), "check")
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

test_that("the shift check scores every unit, and the rounds follow its step", {
  skip_without_julia()
  # One effect a subject on a nonlinear drift: the check finds correcting
  # Laplace would move an estimate past 0.1 se, the units carrying that go to
  # the rounds (on product rules, each one effect wide), and the rounds move
  # the largest coordinate the way the step pointed. On the linear fixture
  # Laplace is exact, and no check runs.
  fits <- .lc_fits("nonlinear")
  corr <- fits$quadrature$laplace$correction
  check <- corr$check
  expect_true(is.list(check))
  expect_identical(check$seeds, 4L)
  expect_length(check$shift, length(fits$quadrature$estimate$raw))
  expect_true(check$status %in% c("flagged", "unresolved"))
  expect_gt(check$carrying, 0L)
  expect_true(isTRUE(corr$applied))
  w <- which.max(abs(check$shift))
  expect_identical(sign(unname(check$shift[w])), sign(unname(corr$delta_se[w])))
  expect_null(.lc_fits("linear")$quadrature$laplace$correction$check)
})

test_that("a corrected fit reports how far its objective is from the exact marginal", {
  skip_without_julia()
  fits <- .lc_fits("linear")
  # Linear Gaussian: Laplace is exact, the screen flags nothing, and the error
  # is zero with nothing to say in print.
  err <- fits$quadrature$laplace$correction$error
  expect_true(is.list(err))
  expect_identical(err$nats, 0)
  expect_identical(err$se, 0)
  expect_identical(err$unscored, 0L)
  expect_false(any(grepl("integrated to within", capture.output(print(fits$quadrature)))))
  # The summary is a sum of its parts, and se is the certification's arithmetic.
  info <- list(error = 0.03, error_draws = 0.01, error_nodes = 0.005,
    error_residual = 0.015, centre_error = c(NaN, 0.01, 0.005))
  e <- .ctLaplaceIntegrationError(info, unscored = 2L)
  expect_equal(e$nats, 0.03)
  expect_equal(e$se, sqrt(0.06))
  expect_equal(e$draws + e$nodes + e$residual, e$nats)
  expect_identical(e$unscored, 2L)
  expect_true(is.na(.ctLaplaceIntegrationError(list(error = NA_real_))$se))
})

test_that("draws' defaults reach the engine, and a flag names both corrections", {
  # No julia needed.
  control <- .ctLaplaceContinueDefaults
  expect_identical(control$draws_se, 0.05)
  expect_identical(control$draws_max, 16384L)
  expect_identical(control$draws_df, 4)
  expect_identical(.ctLaplaceContinuationArgs(control)$draws_df, 4)
  # Units on draws make both routes slow, and the message says so without
  # dropping either.
  plain <- .ctLaplaceShiftMessage(0.63, "popsd_drift", 27L, 40L)
  wide <- .ctLaplaceShiftMessage(2.3, "popsd_mm_Y5", 97L, 100L, drawn = 97L)
  for (said in c(plain, wide)) {
    expect_match(said, "likelihood = 'quadrature'", fixed = TRUE)
    expect_match(said, "intoverpop = 'none'", fixed = TRUE)
  }
  expect_no_match(plain, "hours", fixed = TRUE)
  expect_match(wide, "Both may take hours here, where 97 of those units",
    fixed = TRUE)
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
