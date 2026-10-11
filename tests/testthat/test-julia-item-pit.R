# Probability integral transforms of non-Gaussian items (item_pit.jl), which
# the adequacy families use in place of their errors.
#
# Checked against what they must be: a censored value's bounds are the
# normal's closed form; a binary item's are one minus its success probability
# and one; and under the model that generated the data, the transforms through
# the normal quantile are standard normal and uncorrelated with the same
# item's at the previous row -- on the filter, and on the Laplace route, where
# they come from the exact predictive.

skip_without_julia()

.pit_model <- function(fixed) suppressMessages(ctModel(type = 'ct',
  manifestNames = paste0('Y', 1:4), latentNames = 'eta',
  LAMBDA = if (fixed) matrix(c(1, 1.2, 1, 0.8), 4) else matrix(c(1, 'l2', 'l3', 'l4'), 4),
  MANIFESTMEANS = if (fixed) matrix(c(0.5, 0.2, 0, 0.3), 4) else
    matrix(c('mm1', 'mm2', 0, 'mm4'), 4),
  MANIFESTVAR = if (fixed) diag(c(0.5, 0, 0, 0.4)) else diag(c('mv1', 0, 0, 'mv4')),
  DRIFT = if (fixed) matrix(-0.5) else 'auto', DIFFUSION = if (fixed) matrix(0.8) else 'auto',
  CINT = matrix(0), T0MEANS = if (fixed) matrix(0) else 'auto',
  T0VAR = if (fixed) matrix(1) else 'auto',
  manifesttype = c(0L, 1L, 2L, 4L), ncategories = c(0L, 0L, 4L, 0L),
  censormin = c(-Inf, -Inf, -Inf, 0), silent = TRUE))

.pit_fits <- function() fit_cached('item pit fits', {
  set.seed(5)
  d <- data.frame(suppressMessages(ctGenerate(.pit_model(TRUE), n = 120, Tpoints = 30,
    dtmean = 1, logdtsd = 0.5, backend = 'julia')))
  m2 <- .pit_model(FALSE)
  m2$pars$indvarying <- m2$pars$param %in% 'mm1'
  list(data = d,
    filter = suppressMessages(ctFit(d, .pit_model(FALSE), backend = 'julia', cores = 2,
      verbose = 0)),
    laplace = suppressMessages(ctFit(d, m2, backend = 'julia', intoverpop = 'laplace',
      cores = 2, verbose = 0)))
})

test_that("a censored value's bounds are the normal's, a binary item's its probabilities", {
  module <- ctsem:::.ctJuliaModule(.pit_fits()$filter$model_spec$project)
  cdf <- function(ebar, s2, y, extras, kind) as.numeric(ctsem:::.ctJuliaEval(sprintf(
    paste0('let r = ContinuousTimeSEM._binary_rule(); collect(ContinuousTimeSEM._ctsem_item_cdf(',
      '%.17g, %.17g, %.17g, Float64[%s], %d, r[1], r[2])) end'), ebar, s2, y,
    paste(extras, collapse = ', '), kind)))
  # Censored below at 0, sd 0.4: interior values a density, the floor a mass.
  V <- 0.7 + 0.4^2
  expect_equal(cdf(0.3, 0.7, 1.1, c(0, Inf, 0.4), 4L), rep(pnorm((1.1 - 0.3) / sqrt(V)), 2),
    tolerance = 1e-10)
  expect_equal(cdf(0.3, 0.7, 0, c(0, Inf, 0.4), 4L), c(0, pnorm((0 - 0.3) / sqrt(V))),
    tolerance = 1e-10)
  # Binary: P(0) at y = 0, then from it to one at y = 1.
  b0 <- cdf(0.4, 0.5, 0, numeric(), 1L); b1 <- cdf(0.4, 0.5, 1, numeric(), 1L)
  expect_equal(b0[1], 0)
  expect_equal(b0[2], b1[1])
  expect_equal(b1[2], 1)
  # The probability of a zero, against the logistic-normal integral.
  p0 <- integrate(function(e) (1 - plogis(e)) * dnorm(e, 0.4, sqrt(0.5)), -Inf, Inf)$value
  expect_equal(b0[2], p0, tolerance = 1e-6)
  # Ordinal: increasing bounds that tile [0, 1] over the categories.
  ord <- t(vapply(1:4, function(k) cdf(-0.2, 0.6, k, c(-1, 0.4, 1.5), 2L), numeric(2)))
  expect_equal(ord[1, 1], 0)
  expect_equal(ord[-1, 1], ord[-4, 2], tolerance = 1e-12)
  expect_equal(ord[4, 2], 1, tolerance = 1e-12)
})

test_that("non-Gaussian items' transforms are standard normal under the model", {
  fits <- .pit_fits()
  for (route in c('filter', 'laplace')) {
    err <- ctsem:::.ctAdequacyErrors(fits[[route]])
    expect_false(err$approximate)
    for (j in 2:4) {
      e <- err$E[, j]; e <- e[is.finite(e)]
      expect_gt(length(e), 3000)
      # 3600 rows: the mean within 0.07 and the sd within 0.07 of one are some
      # four standard errors; estimation moves them less than that.
      expect_lt(abs(mean(e)), 0.07, label = paste(route, j, 'mean'))
      expect_lt(abs(stats::sd(e) - 1), 0.07, label = paste(route, j, 'sd'))
      expect_gt(suppressWarnings(stats::ks.test(e, 'pnorm')$p.value), 0.001)
    }
    # Each item's transforms uncorrelated with its own at the subject's
    # previous row.
    same <- which(c(FALSE, err$unit[-1] == err$unit[-length(err$unit)]))
    for (j in 2:4) expect_lt(abs(stats::cor(err$E[same - 1, j], err$E[same, j],
      use = 'complete')), 0.06)
  }
})
