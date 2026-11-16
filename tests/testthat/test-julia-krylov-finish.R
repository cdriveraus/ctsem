# The finish's Krylov model of the Hessian (`CTSEMKrylov`, optimiser.jl): the
# units' score outer product corrected by Hessian-vector products. One Laplace
# fit, made with the model named (`newton = 'krylov'`; by default it is used
# only where the exact Hessian costs 200 gradients or more), and checks against
# the exact Hessian and the chord finish on the same objective.

.krylov_fit <- function() fit_cached("krylov_fit", {
  m <- ctModel(type = "ct", manifestNames = c("Y1", "Y2"),
    latentNames = c("eta1", "eta2"), LAMBDA = diag(2), silent = TRUE)
  m$pars$indvarying <- (m$pars$matrix %in% c("CINT", "T0MEANS") |
    (m$pars$matrix == "DRIFT" & m$pars$row == m$pars$col)) & is.na(m$pars$value)
  set.seed(1)
  suppressWarnings(suppressMessages(ctFit(ctstantestdat, m, backend = "julia",
    intoverpop = "laplace", priors = TRUE, cores = 2,
    optimcontrol = list(newton = "krylov"))))
})

test_that("a fit stepping on the Krylov model ends certified on the exact Hessian", {
  skip_without_julia()
  fit <- .krylov_fit()
  expect_true(fit$optim$converged)
  expect_gt(fit$optim$newton_krylov_models, 0L)
  expect_gt(fit$optim$newton_krylov_products, 0L)
  # The final Hessian is the exact one: within the reuse distance of the
  # estimate, as the certification requires.
  expect_lte(fit$optim$hessian_distance, 0.1)
})

test_that("the Krylov model is the Hessian once its space is the whole space", {
  skip_without_julia()
  fit <- .krylov_fit()
  est <- as.numeric(fit$estimate$raw)
  npar <- length(est)
  model <- ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_krylov_hessian",
    ctsem:::.ctJuliaObjective(fit), est, products = npar, rtol = -Inf)
  Hk <- matrix(as.numeric(model), npar, npar)
  # Forward differences of the gradient along the Lanczos vectors against
  # central differences by coordinate, at the same point: the fit's own
  # Hessian may be from up to a tenth of a standard error away
  # (`hessianreuse`), which moves it by a per cent here.
  exact <- ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_laplace_hessian",
    ctsem:::.ctJuliaObjective(fit), est)
  expect_equal(Hk, matrix(as.numeric(exact), npar, npar), tolerance = 1e-3)
})

test_that("the finish on the Krylov model reaches the chord's maximum for no more Hessians", {
  skip_without_julia()
  fit <- .krylov_fit()
  est <- as.numeric(fit$estimate$raw)
  npar <- length(est)
  V <- solve(-(fit$uncertainty$hessian + t(fit$uncertainty$hessian)) / 2)
  root <- chol((V + t(V)) / 2)
  ctsem:::.ctJuliaEval(paste(
    "function _test_krylov_finish(o, x0, curvature)",
    "  fg! = ContinuousTimeSEM._ctsem_trial_closure(o, :adjoint)",
    "  G = zeros(length(x0)); f = fg!(0.0, G, collect(Float64, x0))",
    "  out = ContinuousTimeSEM._ctsem_newton_finish(o, collect(Float64, x0), f, G,",
    "    fg!; tol=1e-8, certify_tol=1e-6, curvature=Symbol(curvature),",
    "    krylov_min_cost=0.0, probe=false)",
    "  (value=-out.f, hessians=out.full_hessians, products=out.krylov_products)",
    "end", sep = "\n"))
  obj <- ctsem:::.ctJuliaObjective(fit)
  # Two standard errors out in every direction at once: far enough that the
  # first Hessian's walk does not close the gap alone. From several starts,
  # since on this fixture both finishes take two Hessians nearly always and
  # either sometimes a third: from sixteen starts, at two fits of it, the
  # chord took 34 and the model 33 (dev3, 2026-11-01), and which one a single
  # start gives turned on where the fit, at two cores, had stopped.
  #
  # Each ends at a maximum no lower than the fit's, not at the other's value:
  # the objective has a second maximum 0.22 nats above the fit's, and from
  # one of these starts the chord's walk ended there where the model's ended
  # at the fit's (dev4, 2026-11-16; the start moves with the fit's last
  # digits).
  top <- ctsem:::.ctJuliaGet(ctsem:::.ctJuliaCall("Main._test_krylov_finish",
    obj, est, "chord"))$value
  for (seed in 4:7) {
    set.seed(seed)
    start <- est + 2 * as.numeric(t(root) %*% rnorm(npar))
    out <- lapply(c(chord = "chord", krylov = "krylov"), function(curvature)
      ctsem:::.ctJuliaGet(ctsem:::.ctJuliaCall("Main._test_krylov_finish", obj,
        start, curvature)))
    expect_gt(out$krylov$value, top - 1e-6)
    expect_gt(out$chord$value, top - 1e-6)
    expect_gt(out$krylov$products, 0)
    expect_lte(out$krylov$hessians, out$chord$hessians + 1)
  }
})

test_that("the Krylov model's walk gives way to the exact Hessian within that Hessian's reach", {
  skip_without_julia()
  fit <- .krylov_fit()
  est <- as.numeric(fit$estimate$raw)
  npar <- length(est)
  hessian <- (fit$uncertainty$hessian + t(fit$uncertainty$hessian)) / 2
  root <- chol(solve(-hessian))
  ctsem:::.ctJuliaEval(paste(
    "function _test_krylov_reach(o, x0)",
    "  fg! = ContinuousTimeSEM._ctsem_trial_closure(o, :adjoint)",
    "  x = collect(Float64, x0)",
    "  G = zeros(length(x)); f = fg!(0.0, G, x)",
    "  out = ContinuousTimeSEM._ctsem_newton_finish(o, x, f, G, fg!;",
    "    tol=1e-8, certify_tol=1e-6, curvature=:krylov, krylov_min_cost=0.0,",
    "    probe=false)",
    "  (start=-f, value=-out.f, hessians=out.full_hessians,",
    "    models=out.krylov_models, steps=out.steps, at=out.hessian_at,",
    "    x=out.x, distance=out.distance)",
    "end", sep = "\n"))
  obj <- ctsem:::.ctJuliaObjective(fit)
  call <- function(x0) ctsem:::.ctJuliaGet(ctsem:::.ctJuliaCall(
    "Main._test_krylov_reach", obj, x0))
  top <- call(est)$value
  # A start with about 0.001 nats to gain: under the reach of the exact
  # Hessian, a tenth of a standard error or 0.005 nats, and far over the
  # certification's tolerance. The model says the point is that near, so the
  # exact Hessian is formed where the walk starts and its steps close the
  # rest. Held to `tol`, the model's steps closed it first and the Hessian was
  # formed where they ended.
  for (seed in 4:6) {
    set.seed(seed)
    z <- rnorm(npar)
    start <- est + sqrt(2 * 0.001) * as.numeric(t(root) %*% (z / sqrt(sum(z^2))))
    out <- call(start)
    expect_gt(out$models, 0)
    expect_equal(out$hessians, 1)
    expect_equal(as.numeric(out$at), start)
    expect_gte(out$steps, 1)
    expect_gt(out$value - out$start, 5e-4)
    expect_lt(top - out$value, 1e-5)
    expect_lte(out$distance, 0.1)
  }
})

test_that("a finish that is not certified forms no exact Hessian in the Krylov model's place", {
  skip_without_julia()
  fit <- .krylov_fit()
  est <- as.numeric(fit$estimate$raw)
  npar <- length(est)
  hessian <- (fit$uncertainty$hessian + t(fit$uncertainty$hessian)) / 2
  root <- chol(solve(-hessian))
  ctsem:::.ctJuliaEval(paste(
    "function _test_krylov_uncertified(o, x0, final, tol, seeded, hessian0, at0)",
    "  fg! = ContinuousTimeSEM._ctsem_trial_closure(o, :adjoint)",
    "  x = collect(Float64, x0)",
    "  G = zeros(length(x)); f = fg!(0.0, G, x)",
    "  out = ContinuousTimeSEM._ctsem_newton_finish(o, x, f, G, fg!;",
    "    tol=Float64(tol), certify_tol=min(Float64(tol), 1e-6),",
    "    curvature=:krylov, krylov_min_cost=0.0, probe=false,",
    "    final=Symbol(final),",
    "    hessian0=seeded ? Matrix{Float64}(hessian0) : nothing,",
    "    hessian0_at=seeded ? collect(Float64, at0) : nothing)",
    "  (x=out.x, hessians=out.full_hessians, models=out.krylov_models,",
    "    walked=out.walked !== nothing, final=out.hessian !== nothing)",
    "end", sep = "\n"))
  obj <- ctsem:::.ctJuliaObjective(fit)
  call <- function(x0, final, tol, at0 = NULL)
    ctsem:::.ctJuliaGet(ctsem:::.ctJuliaCall("Main._test_krylov_uncertified",
      obj, x0, final, tol, !is.null(at0), hessian,
      if (is.null(at0)) x0 else at0))
  # The first stage of a fit with smoothing parameters (`final = :dear`):
  # however the model's walk ends, it is what the finish returns. Where the
  # model could not be extended, having cost what the exact Hessian would,
  # the exact one was formed and walked on.
  for (seed in 4:7) {
    set.seed(seed)
    start <- est + 2 * as.numeric(t(root) %*% rnorm(npar))
    out <- call(start, "dear", 1e-8)
    expect_gt(out$models, 0)
    expect_equal(out$hessians, 0)
    expect_true(out$walked)
    expect_false(out$final)
  }
  # A pass (`final = :walked`) handed a matrix from elsewhere that cannot make
  # a step: at the maximum to what the value resolves, no step can gain, so
  # the handed matrix's is refused. The model takes over, and its step is
  # refused in turn, which ends the walk. The exact Hessian was formed there.
  top <- call(est, "walked", 1e-13, est)$x
  set.seed(8)
  away <- top + 2 * as.numeric(t(root) %*% rnorm(npar))
  out <- call(top, "walked", 0, away)
  expect_gt(out$models, 0)
  expect_equal(out$hessians, 0)
  expect_true(out$walked)
})
