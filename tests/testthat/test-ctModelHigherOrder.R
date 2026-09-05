# `ctModelHigherOrder` was the one export named in no test file at all
# (review J15, section 9). It is pure model-object surgery, so this needs no
# fit: what it must get right is which cell of which matrix each parameter ends
# up in, and that is exactly what a silent rewrite would get wrong.

.ho_base <- function() suppressMessages(ctModel(LAMBDA = diag(1, 2),
  DRIFT = 0, MANIFESTMEANS = 0, type = 'omx', Tpoints = 4))

test_that("raising the order adds a derivative latent per index", {
  om <- .ho_base()
  ho <- ctModelHigherOrder(om, 1:2)

  expect_equal(ho$n.latent, 4)
  expect_equal(ho$latentNames, c('eta1', 'eta2', 'deta1', 'deta2'))
  expect_equal(dim(ho$DRIFT), c(4L, 4L))
  expect_equal(dim(ho$LAMBDA), c(2L, 4L))
  # The manifests still load on the original latents only.
  expect_equal(unname(ho$LAMBDA[, 3:4]), matrix(0, 2, 2))
})

test_that("the autoregressive effect moves to the derivative", {
  ho <- ctModelHigherOrder(.ho_base(), 1:2)

  # Order-1 self effects are gone, and each level is driven by its derivative
  # with a fixed unit coefficient -- that is what makes it a derivative.
  expect_equal(unname(diag(ho$DRIFT)), rep('0', 4))
  expect_equal(unname(ho$DRIFT[1, 3]), '1')
  expect_equal(unname(ho$DRIFT[2, 4]), '1')
  # The estimated effect is of the level on its own derivative, named for the
  # pair and bounded away from explosive by default.
  expect_equal(unname(ho$DRIFT[3, 1]),
    'drift_deta1_eta1|-log1p(exp(-param*2))-1e-6')
  expect_equal(unname(ho$DRIFT[4, 2]),
    'drift_deta2_eta2|-log1p(exp(-param*2))-1e-6')
  # No cross coupling unless it is asked for.
  expect_equal(unname(ho$DRIFT[3, 2]), '0')
  expect_equal(unname(ho$DRIFT[4, 1]), '0')
})

test_that("explosive = TRUE drops the equilibrium bound", {
  ho <- ctModelHigherOrder(.ho_base(), 1:2, explosive = TRUE)
  expect_equal(unname(ho$DRIFT[3, 1]), 'drift_deta1_eta1|param*2-1')
  expect_equal(unname(ho$DRIFT[4, 2]), 'drift_deta2_eta2|param*2-1')
})

test_that("diffusion follows the process it belongs to", {
  moved <- ctModelHigherOrder(.ho_base(), 1:2)
  # With diffusion = TRUE (the default) the process noise is on the
  # derivatives and the levels are noise free.
  expect_equal(unname(moved$DIFFUSION[1:2, ]), matrix('0', 2, 4))
  expect_equal(unname(moved$DIFFUSION[3, 3]), 'diff_eta1')
  expect_equal(unname(moved$DIFFUSION[4, 3]), 'diff_eta2_eta1')
  expect_equal(unname(moved$DIFFUSION[4, 4]), 'diff_eta2')

  kept <- ctModelHigherOrder(.ho_base(), 1:2, diffusion = FALSE)
  expect_equal(unname(kept$DIFFUSION[1, 1]), 'diff_eta1')
  expect_equal(unname(kept$DIFFUSION[3:4, 3:4]), matrix('0', 2, 2))
})

test_that("initial states and their covariance cover the new latents", {
  ho <- ctModelHigherOrder(.ho_base(), 1:2)
  expect_equal(unname(ho$T0MEANS[, 1]),
    c('T0m_eta1', 'T0m_eta2', 'T0mean_d_eta1', 'T0mean_d_eta2'))
  # T0VAR is filled on and below the diagonal for every new latent against
  # every latent, so the derivatives are free to covary with the levels.
  expect_equal(unname(ho$T0VAR[3, 1]), 'T0var_deta1_eta1')
  expect_equal(unname(ho$T0VAR[4, 2]), 'T0var_deta2_eta2')
  expect_equal(unname(ho$T0VAR[4, 4]), 'T0var_deta2_deta2')
  expect_equal(unname(ho$T0VAR[1, 2]), '0')
})

test_that("the raised model is still a model ctsem will take", {
  ho <- ctModelHigherOrder(.ho_base(), 1:2)
  sm <- suppressMessages(ctModelConvertOMX(ho))
  expect_s3_class(sm, 'ctStanModel')
  pars <- stats::na.omit(sm$pars$param)
  expect_true(all(c('drift_deta1_eta1', 'drift_deta2_eta2',
    'T0mean_d_eta1', 'T0var_deta2_deta2') %in% pars))
})
