# optimcontrol$saem: SAEM before the quasi-Newton optimiser on the Laplace
# route (`ctsem_saem` in the engine).
#
# The engine's own suite checks SAEM against closed forms -- its draws against
# the exact conditional of the random effects, its gradient against the exact
# marginal one, its fixed point against the exact posterior mode
# (test_saem.jl). What is left for here is the R side: that the switch reaches
# the engine on the first stage, that the fit ends where a Laplace fit of the
# same model ends (the averaged point is polished on the Laplace objective and
# certified as usual), that the phase is recorded, that cores > 1 and
# set.seed() behave, and that every way of asking where it cannot apply is
# refused by name.
#
# One SAEM fit of the shared Laplace fixture, cached; the reproducibility check
# pays for a second.

.saem_fit <- function(cores = 1L, seed = 4L) {
  fit_cached(paste("saem_fixture", cores, seed), {
    set.seed(seed)
    suppressWarnings(suppressMessages(ctsem::ctFit(laplace_fixture_data(),
      laplace_fixture_model(), backend = "julia", cores = cores,
      intoverpop = "laplace", priors = TRUE,
      optimcontrol = list(finishsamples = 20, saem = 1000))))
  })
}

test_that("a saem fit ends at the Laplace fit's optimum and records its phase", {
  skip_without_julia()
  fit <- .saem_fit()
  ref <- laplace_fixture()
  # The averaged point is handed to the same optimiser and certification, so
  # where both end at a maximum they end at the same one.
  expect_equal(as.numeric(fit$estimate$loglik), as.numeric(ref$estimate$loglik),
    tolerance = 1e-4)
  expect_equal(as.numeric(fit$estimate$raw), as.numeric(ref$estimate$raw),
    tolerance = 1e-3)
  op <- fit$optim
  expect_gt(op$saem_iterations, 0L)
  expect_lte(op$saem_iterations, 1000L)
  expect_type(op$saem_settled, "logical")
  expect_gte(op$saem_chains, 1L)
  expect_true(is.finite(op$saem_trend))
  expect_true(op$saem_acceptance > 0.05 && op$saem_acceptance < 0.9)
  tr <- op$saem_trace
  expect_s3_class(tr, "data.frame")
  expect_equal(nrow(tr), op$saem_iterations)
  expect_true(all(c("iteration", "logpost_complete", "gradient_norm",
    "step", "acceptance", "trend") %in% names(tr)))
  expect_true(all(is.finite(tr$logpost_complete)))
  # It stopped on its trend rule: no parameter still travelling.
  expect_true(op$saem_settled)
  expect_lte(op$saem_trend, 18 / 16)
  # And a fit that did not ask ran none.
  expect_identical(ref$optim$saem_iterations, 0L)
  expect_null(ref$optim$saem_trace)
})

test_that("saem runs at cores > 1 and set.seed() reproduces it", {
  skip_without_julia()
  a <- .saem_fit(cores = 2L, seed = 5L)
  set.seed(5L)
  b <- suppressWarnings(suppressMessages(ctsem::ctFit(laplace_fixture_data(),
    laplace_fixture_model(), backend = "julia", cores = 2L,
    intoverpop = "laplace", priors = TRUE,
    optimcontrol = list(finishsamples = 20, saem = 1000))))
  expect_identical(a$optim$saem_trace, b$optim$saem_trace)
  expect_equal(as.numeric(a$estimate$raw), as.numeric(b$estimate$raw))
  # And at two cores it ends where one core ends.
  expect_equal(as.numeric(a$estimate$loglik), as.numeric(.saem_fit()$estimate$loglik),
    tolerance = 1e-4)
})

test_that("saem is refused by name where there is nothing for it to sample", {
  dat <- laplace_fixture_data()
  model <- laplace_fixture_model()
  expect_error(ctsem::ctFit(dat, model, backend = "stan", fit = FALSE,
    optimcontrol = list(saem = TRUE)), "optimcontrol\\$saem")
  # FALSE describes what stan does, and is accepted.
  expect_type(suppressWarnings(suppressMessages(ctsem::ctFit(dat, model,
    backend = "stan", fit = FALSE, optimcontrol = list(saem = FALSE)))), "list")
  skip_without_julia()
  expect_error(ctsem::ctFit(dat, model, backend = "julia", fit = FALSE,
    intoverpop = "augmented", optimcontrol = list(saem = TRUE)),
    "optimcontrol\\$saem samples the random effects")
  expect_error(ctsem::ctFit(dat, model, backend = "julia", fit = FALSE,
    intoverpop = "laplace", optimcontrol = list(saem = -3)),
    "TRUE, FALSE or a positive number")
  expect_error(ctsem::ctFit(dat, model, backend = "julia", fit = FALSE,
    intoverpop = "laplace", optimcontrol = list(saem = c(10, 20))),
    "TRUE, FALSE or a positive number")
})
