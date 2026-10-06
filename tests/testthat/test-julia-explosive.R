# Explosive dynamics over a long gap: the engine counts such a subject's filter
# passes (`_CTSEM_EXPLOSIVE_GROWTH`), takes its gradient by forward mode when
# `optimcontrol$explosive_forward` asks (`_CTSEM_ADJOINT_GROWTH`), and the fit
# warns only when the estimate itself, or the posterior, has such a subject.

skip_without_julia()

# A fixed drift of +0.4 and a 20-unit gap: every subject's prediction grows by
# about e^8 over it, so every subject is explosive at any estimate. Simulated
# here rather than by ctGenerate, whose draw stream moves under unrelated
# commits.
.explosive_fit <- function(drift = .4, optimcontrol = list()) {
  set.seed(11)
  times <- c(0, 1, 2, 3, 23)
  d <- do.call(rbind, lapply(1:6, function(i) {
    x <- numeric(length(times))
    x[1] <- stats::rnorm(1)
    for (k in 2:length(times)) {
      x[k] <- exp(.4 * (times[k] - times[k - 1])) * x[k - 1] + stats::rnorm(1, 0, .5)
    }
    data.frame(id = i + 100, time = times,
      Y1 = x + stats::rnorm(length(times), 0, .5))
  }))
  m <- ctModel(type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), DRIFT = matrix(drift), DIFFUSION = matrix(.5),
    MANIFESTVAR = matrix("mvar"), MANIFESTMEANS = matrix("mm"),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0), silent = TRUE)
  m$pars$indvarying <- FALSE
  warnings <- character()
  fit <- withCallingHandlers(
    suppressMessages(ctFit(d, m, backend = "julia", cores = 1,
      optimcontrol = optimcontrol)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(fit = fit, warnings = warnings)
}

test_that("explosive dynamics at the estimate are counted, named and warned about", {
  run <- .explosive_fit()
  fit <- run$fit
  expect_gt(fit$optim$explosive_passes, 0)
  # The fallback is opt-in, so nothing took forward mode.
  expect_equal(fit$optim$forward_gradients, 0)
  expect_equal(sort(as.numeric(fit$optim$explosive_subjects)), 101:106)
  expect_true(any(grepl("explosive dynamics", run$warnings)))

  # The posterior check, on two copies of the estimate as draws.
  draws <- rbind(fit$estimate$raw, fit$estimate$raw)
  expect_warning(ex <- ctsem:::.ctBackendExplosiveDraws(fit, draws, 1L),
    "In 2 of 2 posterior draws")
  expect_equal(ex$draws[["explosive"]], 2L)
  expect_equal(sort(as.numeric(ex$subjects)), 101:106)
})

test_that("a stable drift over the same gap takes no forward gradient and warns nothing", {
  run <- .explosive_fit(drift = -.4)
  expect_equal(run$fit$optim$explosive_passes, 0)
  expect_equal(run$fit$optim$forward_gradients, 0)
  expect_null(run$fit$optim$explosive_subjects)
  expect_false(any(grepl("explosive", run$warnings)))
})

test_that("explosive_forward takes such subjects' gradients by forward mode, for the fit only", {
  run <- .explosive_fit(optimcontrol = list(explosive_forward = TRUE))
  expect_gt(run$fit$optim$forward_gradients, 0)
  expect_equal(sort(as.numeric(run$fit$optim$explosive_subjects)), 101:106)
  # Put back afterwards: a session setting must not outlive the fit.
  expect_equal(as.numeric(ctsem:::.ctJuliaEval(
    "ContinuousTimeSEM._CTSEM_ADJOINT_GROWTH[]")), Inf)
})
