# One accuracy in posterior standard deviations (`optimcontrol$accuracy`,
# `ctFitUncertainty`'s `control$accuracy`): the sampler's stop set from it,
# a method's own account of its error read against it, and the escalation
# to sampling where the account falls short. The rule is tested on hand-built
# records; the one fit checks that the control reaches the stage and the
# sampler gets the stop, with the sampler itself mocked.

test_that("the sampler's stop from an accuracy is the endpoint note's calibration", {
  targets <- ctsem:::.ctAccuracyTargets
  # review/APPROX-endpoint-2026-10-04.md: with the 1.35 margin, eps = 0.5
  # means R-hat <= 1.07 and tail ESS >= ~125; eps = 0.3 means 1.025 and ~350.
  half <- targets(0.5)
  expect_equal(half$minESS, 124)
  expect_equal(half$rhatTarget, 1.066, tolerance = 1e-3)
  third <- targets(0.3)
  expect_equal(third$minESS, 344)
  expect_equal(third$rhatTarget, 1.024, tolerance = 1e-3)
  # Tighter asks for more, and nothing else is accepted.
  expect_gt(targets(0.1)$minESS, targets(0.2)$minESS)
  expect_error(targets(0), "positive")
  expect_error(targets(-1), "positive")
  expect_error(targets("a"), "positive")
})

.accuracy_fit <- function(method = "checked", used = TRUE, meanerror = c(a = 0.1, b = 0.2),
  sdratio = c(a = 1, b = 1.1), stein = TRUE, sampled = FALSE) {
  fit <- list(uncertainty = list(settings = list(method = method),
    details = list(importance_sampling = if (method %in% c("checked", "is")) list(
      used = used, stein = if (stein) list(meanerror = meanerror, sdratio = sdratio)))))
  if (sampled) fit$sample <- list(sampler = "nuts")
  fit
}

test_that("a method's own account of its error is read against the accuracy", {
  short <- ctsem:::.ctAccuracyShortfall
  # Within: every mean within the accuracy and every spread within it.
  ok <- short(.accuracy_fit(), 0.5)
  expect_false(ok$short)
  expect_equal(ok$error, 0.2)
  expect_match(ok$reason, "within")
  # A mean further off than asked.
  off <- short(.accuracy_fit(meanerror = c(a = 0.1, b = 0.7)), 0.5)
  expect_true(off$short)
  expect_equal(off$error, 0.7)
  expect_match(off$reason, "mean")
  # A spread off by more than the accuracy, either way.
  wide <- short(.accuracy_fit(sdratio = c(a = 1, b = 1.8)), 0.5)
  expect_true(wide$short)
  expect_match(wide$reason, "spread")
  narrow <- short(.accuracy_fit(sdratio = c(a = 1, b = 0.6)), 0.5)
  expect_true(narrow$short)
  # The same answer is within a looser accuracy.
  expect_false(short(.accuracy_fit(sdratio = c(a = 1, b = 1.8)), 1)$short)
  # The sampler fell back to the normal: the error is unknown.
  back <- short(.accuracy_fit(used = FALSE), 0.5)
  expect_true(back$short)
  expect_match(back$reason, "normal")
  expect_true(is.na(back$error))
  # No check, no account.
  expect_true(short(.accuracy_fit(stein = FALSE), 0.5)$short)
  # Methods with no account of their own.
  hess <- short(.accuracy_fit(method = "hessian"), 0.5)
  expect_true(hess$short)
  expect_match(hess$reason, "no account")
  expect_true(short(.accuracy_fit(method = "sandwich"), 0.5)$short)
  # A fit already sampled is never short.
  done <- short(.accuracy_fit(sampled = TRUE), 0.5)
  expect_false(done$short)
  expect_equal(done$method, "sample")
})

test_that("the escalation samples to the accuracy only where the pass is short, and records it", {
  seen <- NULL
  testthat::local_mocked_bindings(.ctBackendUncertaintySample = function(fit, control = list(), cores = 1L, verbose = 0, ...) {
    seen <<- control
    fit$sample <- list(sampler = "nuts")
    fit$uncertainty$settings$method <- "sample"
    fit
  }, .package = "ctsem")
  escalate <- ctsem:::.ctUncertaintyEscalate
  # Within: nothing sampled, the reading recorded.
  fit <- escalate(.accuracy_fit(), 0.5)
  expect_null(seen)
  expect_null(fit$sample)
  expect_false(fit$uncertainty$escalation$sampled)
  expect_equal(fit$uncertainty$escalation$accuracy, 0.5)
  expect_equal(fit$uncertainty$escalation$from, "checked")
  # Short: sampled with the stop the accuracy sets, and the fit says so.
  fit <- escalate(.accuracy_fit(meanerror = c(a = 0.9, b = 0.2)), 0.5,
    control = list(chains = 2L))
  expect_equal(seen$minESS, 124)
  expect_equal(seen$rhatTarget, 1.066, tolerance = 1e-3)
  expect_equal(seen$chains, 2L)
  expect_true(fit$uncertainty$escalation$sampled)
  expect_equal(fit$uncertainty$escalation$from, "checked")
  expect_equal(fit$uncertainty$escalation$error, 0.9)
  expect_equal(fit$uncertainty$settings$method, "sample")
  # No accuracy asked: the fit as it was.
  expect_identical(escalate(.accuracy_fit(), NULL), .accuracy_fit())
})

test_that("optimcontrol$accuracy reaches the uncertainty stage, and ctFitUncertainty's control$accuracy the same way", {
  skip_without_julia()
  seen <- NULL
  testthat::local_mocked_bindings(.ctBackendUncertaintySample = function(fit, control = list(), cores = 1L, verbose = 0, ...) {
    seen <<- control
    fit$sample <- list(sampler = "nuts")
    fit$uncertainty$settings$method <- "sample"
    fit
  }, .package = "ctsem")
  set.seed(1)
  dat <- do.call(rbind, lapply(1:10, function(i) {
    eta <- as.numeric(stats::filter(stats::rnorm(6), 0.6, method = "recursive"))
    data.frame(id = i, time = 0:5, Y1 = eta + stats::rnorm(6, 0, .5))
  }))
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    T0MEANS = matrix(0), CINT = matrix(0), MANIFESTMEANS = matrix(0))))
  # 'hessian' has no account of its error, so an accuracy asked for sends the
  # stage to the sampler, with the stop set from it.
  fit <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    cores = 1, verbose = 0, optimcontrol = list(uncertainty = "hessian",
      accuracy = 0.5, sampleControl = list(chains = 2L), carefulfit = FALSE))))
  expect_equal(seen$minESS, 124)
  expect_equal(seen$chains, 2L)
  expect_true(fit$uncertainty$escalation$sampled)
  expect_equal(fit$uncertainty$escalation$from, "hessian")
  # The same through ctFitUncertainty, on the fit before it was sampled.
  seen <- NULL
  plain <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    cores = 1, verbose = 0, optimcontrol = list(uncertainty = "hessian",
      carefulfit = FALSE))))
  expect_null(plain$uncertainty$escalation)
  again <- suppressWarnings(suppressMessages(ctFitUncertainty(plain, "hessian",
    control = list(accuracy = 0.3))))
  expect_equal(seen$minESS, 344)
  expect_true(again$uncertainty$escalation$sampled)
  # Refused where it cannot be read or reached.
  expect_error(ctFitUncertainty(plain, "sample", control = list(accuracy = 0.3)),
    "minESS and rhatTarget")
  expect_error(ctFitUncertainty(plain, "hessian", control = list(accuracy = 0)),
    "positive")
})
