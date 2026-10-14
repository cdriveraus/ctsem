# ctPlotPosterior() on a julia fit. It refused one until 4.0.0; it now draws
# the julia fit's posterior through the transforms summary() reads it with, and
# the full raw prior through the same transforms, so what it plots can be
# checked against summary() and against a transform whose prior is known.

test_that("ctPlotPosterior draws a julia fit's posterior from the draws summary() reads", {
  skip_without_julia()
  model <- suppressMessages(ctModel(type = 'ct', manifestNames = c('Y1', 'Y2'),
    latentNames = c('eta1', 'eta2'), LAMBDA = diag(2)))
  fit <- suppressWarnings(suppressMessages(ctFit(ctstantestdat, model,
    backend = 'julia', cores = 1)))
  set.seed(1)
  pars <- ctsem:::.ctPlotPosteriorDrawsJulia(fit, 4000)
  s <- summary(fit)
  parnames <- vapply(pars, `[[`, '', 'name')
  # One entry per free population parameter, in summary()'s order, so `rows`
  # counts the way the summary lists them.
  expect_equal(parnames, rownames(s$popmeans))

  mm <- pars[[match('mm_Y1', parnames)]]
  expect_equal(mean(mm$quantities[[1]]$post), s$popmeans['mm_Y1', 'mean'],
    tolerance = 1e-3)
  # MANIFESTMEANS is 10 * raw by default, so its prior is normal(0, 10).
  expect_equal(sd(mm$quantities[[1]]$prior), 10, tolerance = 0.05)
  expect_equal(mean(mm$quantities[[1]]$prior), 0, tolerance = 0.5)

  # A varying parameter also gets its population sd, the one summary() reports.
  t0 <- pars[[match('T0m_eta1', parnames)]]
  expect_equal(vapply(t0$quantities, `[[`, '', 'type'), c('Pop. Mean', 'Pop. SD'))
  expect_equal(median(t0$quantities[[2]]$post), s$popsd['T0m_eta1', '50%'],
    tolerance = 0.01)
  expect_true(all(t0$quantities[[2]]$prior > 0))
  # A parameter that does not vary has no sd panel.
  drift <- pars[[match('drift_eta1', parnames)]]
  expect_length(drift$quantities, 1L)

  plots <- ctPlotPosterior(fit, plot = FALSE, priorsamples = 2000)
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
  # Four parameters a page by default, and `npp` is honoured.
  expect_length(plots, ceiling(length(pars) / 4))
  expect_length(ctPlotPosterior(fit, plot = FALSE, priorsamples = 500, npp = 13), 1L)

  expect_error(ctPlotPosterior(fit, cores = 2), 'takes no further arguments')
})

test_that("ctPlotPosterior says what it needs rather than failing inside", {
  model <- suppressMessages(ctModel(type = 'ct', manifestNames = 'Y1',
    latentNames = 'eta1', LAMBDA = matrix(1)))
  # A model on its own failed with "argument 1 is not a vector".
  expect_error(ctPlotPosterior(model), 'needs a fit')
  expect_error(ctPlotPosterior(list(a = 1)), 'a fit from ctFit')

  skip_without_julia()
  model$pars$indvarying <- FALSE
  bare <- suppressWarnings(suppressMessages(ctFit(ctstantestdat, model,
    backend = 'julia', cores = 1, optimcontrol = list(estonly = TRUE))))
  expect_error(ctPlotPosterior(bare), 'ctFitUncertainty')
})
