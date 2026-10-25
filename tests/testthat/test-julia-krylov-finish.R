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
  Hk <- matrix(as.numeric(ctsem:::.ctJuliaGet(model)), npar, npar)
  # Forward differences of the gradient along the Lanczos vectors against the
  # certification's central differences by coordinate.
  expect_equal(Hk, fit$uncertainty$hessian, tolerance = 1e-3)
})

test_that("the finish on the Krylov model reaches the chord's maximum for fewer Hessians", {
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
  set.seed(4)
  # Two standard errors out in every direction at once: far enough that the
  # first Hessian's walk does not close the gap alone.
  start <- est + 2 * as.numeric(t(root) %*% rnorm(npar))
  out <- lapply(c(chord = "chord", krylov = "krylov"), function(curvature)
    ctsem:::.ctJuliaGet(ctsem:::.ctJuliaCall("Main._test_krylov_finish", obj,
      start, curvature)))
  expect_equal(out$krylov$value, out$chord$value, tolerance = 1e-6,
    scale = 1)
  expect_gt(out$krylov$products, 0)
  expect_lte(out$krylov$hessians, out$chord$hessians)
})
