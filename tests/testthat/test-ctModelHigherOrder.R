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

# The type = 'ct' form: the same surgery on the parameter table, which is where
# a ct model keeps everything about a cell. A moved cell keeps its transform
# and its individual variation; a matrix rebuilt from labels would lose both.
.ho_ct <- function() {
  m <- suppressMessages(ctModel(type = 'ct', silent = TRUE,
    manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
    LAMBDA = diag(2)))
  m$pars$indvarying[m$pars$param %in% 'drift_eta2'] <- TRUE
  m
}

test_that("a ct model is raised cell by cell, keeping each moved cell's settings", {
  m <- .ho_ct()
  ho <- ctModelHigherOrder(m, 2)
  expect_s3_class(ho, 'ctStanModel')
  expect_equal(ho$n.latent, 3L)
  expect_equal(ho$latentNames, c('eta1', 'eta2', 'deta2'))
  mats <- ctModelMatrices(ho)
  expect_equal(dim(mats$DRIFT), c(3L, 3L))
  expect_equal(dim(mats$LAMBDA), c(2L, 3L))
  # The level is driven by its derivative, the derivative pulled back by the
  # level, and the level's own autoregression is now the derivative's damping.
  expect_equal(unname(mats$DRIFT[2, 2]), '0')
  expect_equal(unname(mats$DRIFT[2, 3]), '1')
  expect_equal(unname(mats$DRIFT[3, 2]), 'drift_deta2_eta2')
  expect_equal(unname(mats$DRIFT[3, 3]), 'drift_eta2')
  moved <- ho$pars[ho$pars$matrix == 'DRIFT' & ho$pars$row == 3 & ho$pars$col == 3, ]
  orig <- m$pars[m$pars$param %in% 'drift_eta2', ]
  expect_identical(moved$transform, orig$transform)
  expect_true(moved$indvarying)
  # Cross effects stay on the levels unless asked for.
  expect_equal(unname(mats$DRIFT[1, 2]), 'drift_eta1_eta2')
  expect_equal(unname(mats$DRIFT[2, 1]), 'drift_eta2_eta1')
  # Process noise, and its covariance with eta1, move to the derivative.
  expect_equal(unname(mats$DIFFUSION[2, ]), c('0', '0', '0'))
  expect_equal(unname(mats$DIFFUSION[3, 3]), 'diff_eta2')
  expect_equal(unname(mats$DIFFUSION[3, 1]), 'diff_eta2_eta1')
  # New latent: no indicator, its own initial mean and covariances.
  expect_equal(as.numeric(mats$LAMBDA[, 3]), c(0, 0))
  expect_equal(unname(mats$T0MEANS[3, 1]), 'T0mean_deta2')
  expect_equal(unname(mats$T0VAR[3, ]), c('T0var_deta2_eta1', 'T0var_deta2_eta2', 'T0var_deta2'))
})

test_that("a ct model raised at a latent that is not the first puts each cell on the new row", {
  m <- .ho_ct()
  k <- which(m$pars$matrix == 'CINT' & m$pars$row == 2)
  m$pars$param[k] <- 'cint2'; m$pars$value[k] <- NA; m$pars$transform[k] <- 'param'
  ho <- ctModelHigherOrder(m, 2, cint = TRUE, diffusion = FALSE)
  mats <- ctModelMatrices(ho)
  expect_equal(unname(mats$CINT[, 1]), c('0', '0', 'cint2'))
  expect_equal(unname(mats$DIFFUSION[2, 2]), 'diff_eta2')
  expect_equal(unname(mats$DIFFUSION[3, 3]), '0')
  expect_error(ctModelHigherOrder(.ho_ct(), 3), 'indices')
})
