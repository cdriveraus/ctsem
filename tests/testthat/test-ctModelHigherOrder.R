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

# ctModelSlowComponent: the same kind of surgery, adding a slow process that
# drives a latent through its drift (the SNSF pilot's specification).
test_that("a slow process drives its latent through the drift and nothing else", {
  m <- .ho_ct()
  m$tipredDefault <- TRUE
  sm <- ctModelSlowComponent(m, 1)
  expect_equal(sm$latentNames, c('eta1', 'eta2', 'slow_eta1'))
  mats <- ctModelMatrices(sm)
  # Its persistence and its effect on the latent are free; nothing reaches it
  # from the fast states; it has no indicators and no intercept.
  expect_equal(unname(mats$DRIFT[3, 3]), 'drift_slow_eta1')
  expect_equal(unname(mats$DRIFT[1, 3]), 'drift_eta1_slow_eta1')
  expect_equal(unname(mats$DRIFT[3, 1:2]), c('0', '0'))
  expect_equal(unname(mats$DRIFT[2, 3]), '0')
  expect_equal(as.numeric(mats$LAMBDA[, 3]), c(0, 0))
  expect_equal(as.numeric(mats$CINT[3, 1]), 0)
  # Unit innovations, independent of the fast ones.
  expect_equal(unname(mats$DIFFUSION[3, ]), c('0', '0', '1'))
  # It starts from its stationary mean and variance, and from a free
  # covariance with the latent it drives; zero with the other.
  expect_equal(unname(mats$T0VAR[3, ]), c('T0var_slow_eta1_eta1', '0', 'stationary'))
  expect_identical(unname(mats$T0MEANS[3, 1]), 'stationary')
  p <- sm$pars
  slowcells <- p[(p$row == 3 & p$matrix %in% c('DRIFT', 'T0MEANS', 'T0VAR')) |
    (p$col == 3 & p$matrix == 'DRIFT'), ]
  free <- slowcells[!is.na(slowcells$param), ]
  # No individual variation and no predictor effects on any slow cell, even
  # where the model's default would add them.
  expect_false(any(free$indvarying))
  for (e in grep('_effect$', names(p), value = TRUE)) expect_true(all(free[[e]] == 'FALSE'))
  # The persistence takes ctsem's default drift transform.
  expect_identical(free$transform[free$param == 'drift_slow_eta1'],
    .ctModelDefaultFreePar('DRIFT', 3L, 3L, TRUE)$transform)
  # The latents it was not asked for are untouched.
  cells <- function(q) { q <- q[q$matrix == 'DRIFT' & q$row <= 2 & q$col <= 2, ]
    q$param[order(q$row, q$col)] }
  expect_identical(cells(p), cells(m$pars))
})

test_that("slow processes' innovations correlate with each other, added together or later", {
  both <- ctModelMatrices(ctModelSlowComponent(.ho_ct(), 1:2))
  expect_equal(unname(both$DIFFUSION[4, ]), c('0', '0', 'diff_slow_eta2_slow_eta1', '1'))
  later <- ctModelSlowComponent(ctModelSlowComponent(.ho_ct(), 1), 2)
  expect_equal(later$latentNames, c('eta1', 'eta2', 'slow_eta1', 'slow_eta2'))
  expect_equal(unname(ctModelMatrices(later)$DIFFUSION[4, ]),
    c('0', '0', 'diff_slow_eta2_slow_eta1', '1'))
  expect_error(ctModelSlowComponent(ctModelSlowComponent(.ho_ct(), 1), 1), 'already has')
  expect_error(ctModelSlowComponent(.ho_ct(), 3), 'indices')
})

test_that("the search starts a slow process away from the stationary point", {
  m <- .ho_ct()
  # The effect at one; the persistence a tenth of the latent's, through the
  # inverse of the default transform. Without a fit to read the latent's drift
  # from, the pilot's -0.2.
  st <- ctsem:::.ctImproveSlowStart(m, NULL, 1)
  expect_equal(unname(st['drift_eta1_slow_eta1']), 1)
  tf <- function(param) eval(parse(text = .ctModelDefaultFreePar('DRIFT', 3L, 3L, TRUE)$transform),
    list(param = param, log1p_exp = ctsem:::log1p_exp))
  expect_equal(tf(unname(st['drift_slow_eta1'])), -0.2, tolerance = 1e-6)
  fixedfast <- m
  k <- which(fixedfast$pars$matrix == 'DRIFT' & fixedfast$pars$row == 1 & fixedfast$pars$col == 1)
  fixedfast$pars$param[k] <- NA; fixedfast$pars$value[k] <- -3
  expect_equal(tf(unname(ctsem:::.ctImproveSlowStart(fixedfast, NULL, 1)['drift_slow_eta1'])),
    -0.3, tolerance = 1e-6)
  # The evidence: the latent's own updates still positively correlated over
  # some gap; a negative correlation is not evidence of a slow process.
  ad <- list(statistics = data.frame(family = 'latent', statistic = c('eta1~eta1 gap <=1',
    'eta1~eta1 gap >5', 'eta2~eta2 gap >5', 'eta1~eta2 gap >5'),
    z = c(1, 3, -4, 5), p = c(0.3, 0.003, 6e-5, 6e-7), stringsAsFactors = FALSE))
  # The smallest p, times the gap bins it is the smallest of.
  expect_equal(ctsem:::.ctImproveSlowEvidence(ad, 'eta1'), 0.006)
  expect_true(is.na(ctsem:::.ctImproveSlowEvidence(ad, 'eta2')))
  expect_true(is.na(ctsem:::.ctImproveSlowEvidence(NULL, 'eta1')))
})

# ctModelPeriodicComponent: a clock pair that turns at the period, with no
# noise, set from the subject's first time, driving the latents asked for.
test_that("a cycle is a noiseless clock pair driving the latents asked for", {
  m <- .ho_ct()
  cm <- ctModelPeriodicComponent(m, period = 24, indices = 1)
  expect_equal(cm$latentNames, c('eta1', 'eta2', 'cycle_cos', 'cycle_sin'))
  expect_equal(cm$TIpredNames, c('cycle_cos0', 'cycle_sin0'))
  expect_equal(cm$n.TIpred, 2L)
  mats <- ctModelMatrices(cm)
  w <- 2 * pi / 24
  expect_equal(as.numeric(mats$DRIFT[3:4, 3:4]), c(0, w, -w, 0))
  expect_equal(unname(mats$DRIFT[1, 3:4]), c('drift_eta1_cycle_cos', 'drift_eta1_cycle_sin'))
  expect_equal(unname(mats$DRIFT[2, 3:4]), c('0', '0'))
  expect_equal(unname(mats$DRIFT[3:4, 1:2]), matrix('0', 2, 2))
  expect_true(all(as.numeric(mats$DIFFUSION[3:4, ]) == 0))
  expect_equal(as.numeric(mats$LAMBDA[, 3:4]), rep(0, 4))
  # The initial state reads the phase predictors; each subject may deviate
  # from the clock by one shared sd, independently of the other states.
  expect_equal(unname(mats$T0MEANS[3:4, 1]), c('cycle_cos0 * 1', 'cycle_sin0 * 1'))
  expect_equal(unname(mats$T0VAR[3:4, 3:4]), matrix(c('cycle_T0sd', '0', '0', 'cycle_T0sd'), 2))
  expect_equal(unname(mats$T0VAR[3:4, 1:2]), matrix('0', 2, 2))
  # No predictor reaches any parameter: the phases are read only by T0MEANS.
  p <- cm$pars
  expect_true(all(p$cycle_cos0_effect == 'FALSE') && all(p$cycle_sin0_effect == 'FALSE'))
  # Without individual deviations the clock is exact.
  fixedm <- ctModelMatrices(ctModelPeriodicComponent(m, 24, individual = FALSE))
  expect_true(all(as.numeric(fixedm$T0VAR[3:4, 3:4]) == 0))
  expect_error(ctModelPeriodicComponent(cm, 24), 'already has')
  expect_error(ctModelPeriodicComponent(m, -1), 'period')
})

test_that("the phase predictors are the clock at each subject's first row", {
  d <- data.frame(id = rep(1:2, each = 3), time = c(8, 12, 20, 30, 34, 42))
  out <- ctPeriodicPhases(d, period = 24)
  expect_equal(out$cycle_cos0, rep(cos(2 * pi * c(8, 30) / 24), each = 3))
  expect_equal(out$cycle_sin0, rep(sin(2 * pi * c(8, 30) / 24), each = 3))
  expect_error(ctPeriodicPhases(d, 24, id = 'subject'), 'subject')
})

# ctModelUniqueProcess: a persistent process loading on one item alone.
test_that("an item's own process loads on it alone and has its own dynamics", {
  m <- suppressMessages(ctModel(type = 'ct', silent = TRUE,
    manifestNames = c('Y1', 'Y2', 'Y3'), latentNames = 'eta',
    LAMBDA = matrix(c(1, 'l2', 'l3'), 3)))
  um <- ctModelUniqueProcess(m, c('Y2', 'Y3'))
  expect_equal(um$latentNames, c('eta', 'u_Y2', 'u_Y3'))
  mats <- ctModelMatrices(um)
  expect_equal(as.numeric(mats$LAMBDA[, 2]), c(0, 1, 0))
  expect_equal(as.numeric(mats$LAMBDA[, 3]), c(0, 0, 1))
  expect_equal(unname(diag(mats$DRIFT)[2:3]), c('drift_u_Y2', 'drift_u_Y3'))
  expect_equal(unname(mats$DRIFT[1, 2:3]), c('0', '0'))
  expect_equal(unname(mats$DRIFT[2, c(1, 3)]), c('0', '0'))
  expect_equal(unname(diag(mats$DIFFUSION)[2:3]), c('diff_u_Y2', 'diff_u_Y3'))
  expect_equal(unname(mats$DIFFUSION[2, c(1, 3)]), c('0', '0'))
  expect_equal(as.numeric(mats$CINT[2:3, 1]), c(0, 0))
  expect_equal(as.numeric(mats$T0MEANS[2:3, 1]), c(0, 0))
  # Each starts from its stationary distribution: the variance its persistence
  # and innovations imply.
  expect_equal(unname(diag(mats$T0VAR)[2:3]), c('stationary', 'stationary'))
  expect_equal(ctModelMatrices(ctModelUniqueProcess(m, 2))$LAMBDA[, 2], mats$LAMBDA[, 2])
  expect_error(ctModelUniqueProcess(um, 'Y2'), 'already has')
  expect_error(ctModelUniqueProcess(m, 'Y9'), 'manifest')
})

test_that("items can share a process, loading on the first at one and the rest freely", {
  m <- suppressMessages(ctModel(type = 'ct', silent = TRUE,
    manifestNames = c('Y1', 'Y2', 'Y3'), latentNames = 'eta',
    LAMBDA = matrix(c(1, 'l2', 'l3'), 3)))
  sm <- ctModelUniqueProcess(m, c('Y2', 'Y3'), shared = TRUE)
  expect_equal(sm$latentNames, c('eta', 'u_Y2_Y3'))
  mats <- ctModelMatrices(sm)
  expect_equal(unname(mats$LAMBDA[, 2]), c('0', '1', 'lambda_Y3_u_Y2_Y3'))
  expect_equal(unname(mats$DRIFT[2, ]), c('0', 'drift_u_Y2_Y3'))
  expect_equal(unname(mats$DIFFUSION[2, ]), c('0', 'diff_u_Y2_Y3'))
  expect_equal(unname(mats$T0VAR[2, 2]), 'stationary')
  expect_error(ctModelUniqueProcess(m, 'Y2', shared = TRUE), 'two items')
})

# The search reads errors correlated within an occasion as a process the items
# share, proposed at the structural stage from the pairs family; the score
# test's error covariance is freed only when 'errorcov' asks for it.
test_that("correlated errors get a shared process from the checks; errorcov frees the cell", {
  m <- .ho_ct()
  st <- list(fixed = list(cells = data.frame(cell = 'MANIFESTVAR[Y2,Y1]',
    matrix = 'MANIFESTVAR', row = 2L, col = 1L, LM = 30, p = 1e-8, estimate = -0.2,
    se = 0.04, suggested = 'manifestvar_Y2_Y1', transform = 'param', p_adj = 1e-8,
    class = 1L, stringsAsFactors = FALSE)))
  expect_length(ctsem:::.ctImproveProposals(st, m, 0.05, NULL,
    relax = c('heterogeneity', 'fixed', 'unique')), 0L)
  pe <- ctsem:::.ctImproveProposals(st, m, 0.05, NULL, relax = c('fixed', 'errorcov'))
  expect_identical(pe[[1]]$family, 'fixed')
  expect_identical(pe[[1]]$label, 'MANIFESTVAR[Y2,Y1] freed')
  ad <- list(statistics = data.frame(family = c('pairs', 'residuals'),
    statistic = c('Y1:Y2', 'Y1~Y1 gap <=3'), z = c(-6, -1), p = c(2e-9, 0.3),
    stringsAsFactors = FALSE))
  pr <- ctsem:::.ctImproveStructural(m, 'unique', NULL, ad, NULL, 0.05)
  expect_length(pr, 1L)
  expect_identical(pr[[1]]$label, 'Y1 and Y2 given a process they share')
  expect_equal(pr[[1]]$p, 2e-9)
  expect_true('u_Y1_Y2' %in% pr[[1]]$edit(m, NULL, NULL)$model$latentNames)
  expect_length(ctsem:::.ctImproveStructural(m, c('unique', 'errorcov'), NULL, ad, NULL,
    0.05), 0L)
  ad$statistics$p[1] <- 0.06
  expect_length(ctsem:::.ctImproveStructural(m, 'unique', NULL, ad, NULL, 0.05), 0L)
})

test_that("the fixed-cell screen leaves a latent without innovations as written", {
  cm <- ctModelPeriodicComponent(.ho_ct(), 24, indices = 1, individual = FALSE)
  terms <- ctsem:::.ctScoreFixedModel(cm)$terms
  # No drift into the clock, no covariance with it; its effects on the
  # latents it does not yet drive are candidates.
  expect_false(any(terms$matrix == 'DRIFT' & terms$row %in% 3:4))
  expect_false(any(terms$matrix == 'DIFFUSION' & (terms$row %in% 3:4 | terms$col %in% 3:4)))
  expect_true(any(terms$matrix == 'DRIFT' & terms$row == 2 & terms$col == 3))
  expect_equal(ctsem:::.ctImproveNoisyLatents(cm), 1:2)
})

test_that("the search reads a cycle's period and an item's own correlation from the checks", {
  ad <- list(families = data.frame(family = c('periodic', 'residuals'), p = c(1e-4, 0.01)),
    statistics = data.frame(family = c('periodic', 'periodic', 'periodic', 'residuals',
      'residuals'), statistic = c('Y1 period 24', 'Y2 period 12', 'latent eta1 period 8',
      'Y2~Y2 gap <=3', 'Y1~Y1 gap <=3'), z = c(4, -6, 3, 3.5, -2),
      p = c(6e-5, 2e-9, 3e-3, 5e-4, 0.05), stringsAsFactors = FALSE))
  cyc <- ctsem:::.ctImprovePeriodicEvidence(ad)
  expect_equal(cyc$period, 24)
  # The strongest positive period's p times the periods scanned, not the
  # family's p, so that a cycle ranks against item processes on their terms.
  expect_equal(cyc$p, 3 * 6e-5)
  expect_equal(ctsem:::.ctImproveOwnEvidence(ad, 'residuals', 'Y2'), 5e-4)
  expect_true(is.na(ctsem:::.ctImproveOwnEvidence(ad, 'residuals', 'Y1')))
  expect_null(ctsem:::.ctImprovePeriodicEvidence(NULL))
})

# The volatility terms: a latent's innovation sd scaled by exp(v * state),
# the diagonal cell's own parameter moved to PARS with everything it carries.
test_that("an innovation sd is scaled by a state through PARS, other cells untouched", {
  m <- .ho_ct()
  p0 <- m$pars
  cell0 <- p0[p0$matrix == 'DIFFUSION' & p0$row == 2 & p0$col == 2, ]
  v <- ctsem:::.ctAddDiffusionTerms(m, c(2L, 2L), c('va', 'vb'), c(2L, 1L))
  p <- v$pars
  cell <- p[p$matrix == 'DIFFUSION' & p$row == 2 & p$col == 2, ]
  expect_equal(cell$param, paste0(cell0$param, ' * exp(va * eta2 + vb * eta1)'))
  expect_true(is.na(cell$transform))
  moved <- p[p$matrix == 'PARS' & p$param == cell0$param, ]
  expect_equal(nrow(moved), 1L)
  expect_identical(moved$transform, cell0$transform)
  terms <- p[p$matrix == 'PARS' & p$param %in% c('va', 'vb'), ]
  expect_equal(terms$transform, c('param', 'param'))
  expect_false(any(terms$indvarying))
  # Every other DIFFUSION cell as it was.
  others <- function(q) { q <- q[q$matrix == 'DIFFUSION' & !(q$row == 2 & q$col == 2), ]
    q[order(q$row, q$col), c('row', 'col', 'param', 'value', 'transform')] }
  expect_equal(others(p), others(p0), ignore_attr = TRUE)
  # A fixed sd keeps its value as a factor.
  f <- m
  k <- which(f$pars$matrix == 'DIFFUSION' & f$pars$row == 1 & f$pars$col == 1)
  f$pars$param[k] <- NA; f$pars$value[k] <- 0.5
  fv <- ctsem:::.ctAddDiffusionTerms(f, 1L, 'vc', 1L)$pars
  expect_equal(fv$param[fv$matrix == 'DIFFUSION' & fv$row == 1 & fv$col == 1],
    '0.5 * exp(vc * eta1)')
  # The screen tests every latent with innovations against every state.
  vm <- ctsem:::.ctScoreVolatilityModel(m)
  expect_equal(nrow(vm$terms), 4L)
  expect_setequal(paste(vm$terms$i, vm$terms$m), c('1 1', '1 2', '2 1', '2 2'))
})

test_that("a model's cycles are not proposed again, nor given slow processes or a second order", {
  cm <- ctModelPeriodicComponent(.ho_ct(), 24, indices = 1:2, individual = FALSE)
  expect_equal(ctsem:::.ctModelCyclePeriods(cm), 24)
  expect_equal(ctsem:::.ctModelCyclePeriods(.ho_ct()), numeric())
  ad <- list(families = data.frame(family = 'periodic', p = 1e-4),
    statistics = data.frame(family = 'periodic',
      statistic = c('latent eta1 period 24', 'latent eta1 period 12'), z = c(8, 4),
      p = c(1e-15, 3e-5), stringsAsFactors = FALSE))
  # The strongest period the model lacks, not the one it has.
  expect_equal(ctsem:::.ctImprovePeriodicEvidence(ad)$period, 24)
  expect_equal(ctsem:::.ctImprovePeriodicEvidence(ad, 24)$period, 12)
  expect_null(ctsem:::.ctImprovePeriodicEvidence(ad, c(24, 12)))
  pr <- ctsem:::.ctImproveStructural(cm, c('periodic', 'slow', 'order'), NULL, ad, NULL, 0.05)
  labels <- vapply(pr, `[[`, '', 'label')
  expect_true('a cycle of period 12' %in% labels)
  expect_false('a cycle of period 24' %in% labels)
  expect_false(any(grepl('cycle_(cos|sin)', labels)))
  expect_true(all(c('eta1 driven by a slow process', 'eta1 made second order') %in% labels))
})

test_that("a second order's derivative is given no slow process and no second order", {
  # Raising the derivative's order is a third order by another name.
  m2 <- ctModelHigherOrder(.ho_ct(), 1L)
  pr <- ctsem:::.ctImproveStructural(m2, c('slow', 'order'), NULL, NULL, NULL, 0.05)
  labels <- vapply(pr, `[[`, '', 'label')
  derivative <- setdiff(m2$latentNames, .ho_ct()$latentNames)
  expect_length(derivative, 1L)
  expect_false(any(startsWith(labels, paste0(derivative, ' '))))
  expect_false('eta1 made second order' %in% labels)
  expect_true(all(c('eta2 made second order', 'eta2 driven by a slow process') %in% labels))
})

test_that("the state dependence screen writes no terms into a clock's drift", {
  cm <- ctModelPeriodicComponent(.ho_ct(), 24, indices = 1:2, individual = FALSE)
  terms <- ctsem:::.ctScoreStateModel(cm)$terms
  expect_setequal(unique(terms$equation), c('eta1', 'eta2'))
  expect_true(any(grepl('cycle_cos', terms$by)))
})
