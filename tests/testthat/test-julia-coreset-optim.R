# The optimiser's coreset stages on the Laplace route (`_ctsem_coreset_plan`
# in optimiser.jl): stretches of the L-BFGS stage on a weighted few of the
# units, each checked by a full gradient, with the full objective finishing.
# Two fits of one model from the default start under one seed, with them
# (`optimcontrol$coreset = TRUE`; the default, 'auto', leaves them off for data
# this small) and without: the same optimum and Hessian, in fewer full
# gradients.
# Three hundred subjects of the shared Laplace fixture's design, so that its
# seven parameters leave a coreset a small share of the units. Batching is off
# in both: from the default start it does most of the climbing on this model,
# and after it the two fits took the same full gradients.

.coreset_optim_fits <- function() {
  fit_cached("coreset optimiser fits", {
    data <- laplace_fixture_data(nsub = 300L, tp = 6L, seed = 41L)
    fit_one <- function(coreset) {
      set.seed(1)
      suppressWarnings(suppressMessages(ctsem::ctFit(data, laplace_fixture_model(),
        backend = "julia", cores = 1, intoverpop = "laplace", priors = TRUE,
        optimcontrol = c(list(finishsamples = 20, batch = FALSE),
          if (coreset) list(coreset = TRUE)))))
    }
    list(on = fit_one(TRUE), off = fit_one(FALSE))
  })
}

test_that("coreset stages reach the full objective's optimum in fewer full gradients", {
  skip_on_cran()
  skip_if_not(isTRUE(ctsem::ctJuliaStatus()$available))
  fits <- .coreset_optim_fits()
  on <- fits$on
  off <- fits$off
  # What ran: stages when asked for, none by default.
  expect_gt(sum(on$optim$coreset_sizes), 0)
  expect_gt(on$optim$coreset_iterations, 0)
  expect_equal(sum(off$optim$coreset_sizes), 0)
  expect_equal(off$optim$coreset_iterations, 0L)
  # The finish and its Hessian are the full objective's either way.
  expect_true(isTRUE(on$uncertainty$certification$certified))
  expect_true(isTRUE(off$uncertainty$certification$certified))
  expect_equal(on$estimate$logposterior, off$estimate$logposterior, tolerance = 1e-8)
  H <- off$uncertainty$hessian
  se <- sqrt(diag(solve(-H)))
  expect_equal(as.numeric(on$estimate$raw) / se, as.numeric(off$estimate$raw) / se,
    tolerance = 1e-3)
  expect_equal(on$uncertainty$hessian, H, tolerance = 1e-5)
  # And the point of it.
  expect_lt(on$optim$full_gradients, off$optim$full_gradients)
})

test_that("the coreset weights are non-negative least squares", {
  skip_on_cran()
  skip_if_not(isTRUE(ctsem::ctJuliaStatus()$available))
  # Karush-Kuhn-Tucker at the engine's solution: non-negative, the residual's
  # gradient zero where a weight is positive and pointing outward where it is
  # zero.
  set.seed(2)
  A <- matrix(rnorm(80 * 12), 80)
  # Half the columns enter negatively, so their weights must be held at zero.
  b <- drop(A %*% rep(c(1, -1), each = 6)) + rnorm(80, 0, 0.1)
  x <- as.numeric(ctsem:::.ctBackendJuliaValue(ctsem:::.ctJuliaCall(
    "ContinuousTimeSEM._ctsem_nnls", crossprod(A), ctsem:::.ctJuliaVector(drop(crossprod(A, b))))))
  g <- drop(crossprod(A, A %*% x - b))
  expect_true(all(x >= 0))
  expect_true(any(x == 0) && any(x > 0))
  expect_lt(max(abs(g[x > 0])), 1e-6)
  expect_gt(min(g[x == 0]), -1e-6)
})

# Which fits 'auto' runs them on (`.ctJuliaCoresetOn()`): at least 200
# top-level units and 10,000 data cells, read off the specification; TRUE and
# FALSE as given; never off the Laplace route.
test_that("coreset = 'auto' runs the stages only on large data", {
  on <- ctsem:::.ctJuliaCoresetOn
  spec <- function(units, rows, manifests) list(
    laplace = list(levels = list(list(ngroups = 3L * units), list(ngroups = units))),
    manifest_data = matrix(0, manifests, rows))
  big <- spec(400L, 5000L, 3L)
  expect_true(on(list(), big))
  expect_true(on(list(coreset = "auto"), big))
  expect_false(on(list(coreset = FALSE), big))
  # The outermost level's groups are the units: 1200 subjects in 400 groups
  # pass, 600 subjects in 100 groups do not.
  expect_false(on(list(), spec(100L, 5000L, 3L)))
  expect_false(on(list(), spec(800L, 3000L, 1L)))
  expect_true(on(list(coreset = TRUE), spec(40L, 240L, 1L)))
  expect_false(on(list(), list(manifest_data = matrix(0, 3, 5000))))
})
