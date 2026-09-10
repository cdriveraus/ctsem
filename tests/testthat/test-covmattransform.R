# covmattransform reaching the engine, and meaning the same thing on both
# backends.
#
# This setting has a history: it was validated in R and then ignored by the
# julia path, which reads as support -- a 'cholesky' model sat 3.76 log units
# from the same model on stan with no error and no warning. So the assertions
# here are mostly about the setting actually changing the answer, not about the
# answer itself.

if (requireNamespace('testthat', quietly = TRUE)) {

  covmat_model <- function(tf, k = 2L) {
    DR <- matrix(0, k, k); diag(DR) <- paste0('dr', 1:k, '|-log1p_exp(param)')
    DF <- matrix(0, k, k)
    for (i in 1:k) for (j in 1:i) {
      DF[i, j] <- if (i == j) paste0('df', i, '|log1p_exp(param)')
      else paste0('dfc', i, '_', j)
    }
    MV <- matrix(0, k, k); diag(MV) <- paste0('mv', 1:k, '|log1p_exp(param)')
    m <- suppressMessages(ctModel(type = 'ct', n.latent = k, n.manifest = k,
      LAMBDA = diag(k), manifestNames = paste0('Y', 1:k),
      latentNames = paste0('eta', 1:k), DRIFT = DR, DIFFUSION = DF,
      MANIFESTVAR = MV, CINT = matrix(0, k, 1), MANIFESTMEANS = matrix(0, k, 1),
      T0MEANS = matrix(0, k, 1), T0VAR = diag(1, k)))
    m$covmattransform <- tf
    m
  }
  covmat_data <- function(k = 2L, nsub = 12L, nt = 6L) {
    set.seed(4)
    d <- data.frame(id = rep(seq_len(nsub), each = nt),
      time = rep(seq_len(nt) - 1, nsub))
    for (i in 1:k) d[[paste0('Y', i)]] <- rnorm(nsub * nt)
    d
  }

  test_that("the covariance construction code is the same on both backends", {
    expect_equal(ctsem:::.ctCovMatCode(covmat_model('rawcorr')), 0L)
    expect_equal(ctsem:::.ctCovMatCode(covmat_model('z')), 2L)
    # and the stan path derives the same code from the same setting
    for (tf in c('rawcorr', 'z')) {
      m <- covmat_model(tf)
      expect_equal(ctsem:::.ctCovMatCode(m),
        if (identical(tf, 'z')) 2L else 0L)
    }
  })

  test_that("an unsupported covmattransform is refused by name", {
    skip_on_cran()
    expect_error(suppressMessages(ctFit(datalong = covmat_data(),
      model = covmat_model('cholesky'), backend = 'julia', fit = FALSE,
      cores = 1L)), 'cholesky')
  })

  test_that("covmattransform reaches the julia spec", {
    skip_on_cran()
    for (tf in c('rawcorr', 'z')) {
      p <- suppressWarnings(suppressMessages(ctFit(datalong = covmat_data(),
        model = covmat_model(tf), backend = 'julia', fit = FALSE, cores = 1L)))
      expect_equal(p$covmatcode, if (identical(tf, 'z')) 2L else 0L)
    }
  })

  # The one that would have caught the bug this test file was written for: the
  # objective cache keyed on everything except the construction, so a 'z' model
  # built after a 'rawcorr' one was handed the earlier objective and returned
  # the earlier likelihood, identical to eight decimal places.
  test_that("the objective cache distinguishes the constructions", {
    skip_on_cran()
    p0 <- suppressWarnings(suppressMessages(ctFit(datalong = covmat_data(),
      model = covmat_model('rawcorr'), backend = 'julia', fit = FALSE,
      cores = 1L)))
    pz <- suppressWarnings(suppressMessages(ctFit(datalong = covmat_data(),
      model = covmat_model('z'), backend = 'julia', fit = FALSE, cores = 1L)))
    expect_false(identical(ctsem:::.ctJuliaObjectiveKey(p0),
      ctsem:::.ctJuliaObjectiveKey(pz)))
  })

  test_that("covmattransform='z' changes the likelihood and keeps its gradient", {
    skip_on_cran()
    dat <- covmat_data()
    got <- list()
    for (tf in c('rawcorr', 'z')) {
      p <- suppressWarnings(suppressMessages(ctFit(datalong = dat,
        model = covmat_model(tf), backend = 'julia', fit = FALSE, cores = 1L)))
      np <- ctsem:::.ctBackendNpar(p)
      set.seed(7); pars <- rnorm(np, 0, 0.05)
      g <- ctsem:::ctJuliaEvaluate(p, pars = pars, gradient = TRUE)
      got[[tf]] <- g$value
      # the reverse pass has to honour the same setting the forward did, or the
      # gradient describes a different model than the likelihood
      ll <- function(q) ctsem:::ctJuliaEvaluate(p, pars = q,
        gradient = FALSE)$value
      h <- 1e-5
      idx <- unique(round(seq(1, np, length.out = 4)))
      fd <- vapply(idx, function(i) {
        a <- pars; a[i] <- a[i] + h; b <- pars; b[i] <- b[i] - h
        (ll(a) - ll(b)) / (2 * h) }, numeric(1))
      expect_equal(g$gradient[idx], fd, tolerance = 1e-4)
    }
    expect_false(isTRUE(all.equal(got$rawcorr, got$z)))
  })

  # A fit or model saved before the (-1, 1) map moved into the covariance
  # construction still carries it in its own parameter table, so it would be
  # applied twice: correlations come back shrunk towards zero with nothing
  # raised. That is how it went unnoticed on a shipped fixture here, so the
  # detection is asserted rather than trusted.
  test_that('a covariance off-diagonal carrying the old transform is reported', {
    skip_on_cran()
    m <- covmat_model('rawcorr')
    off <- m$pars$matrix %in% c('DIFFUSION', 'MANIFESTVAR', 'T0VAR') &
      m$pars$row != m$pars$col & is.na(m$pars$value)
    expect_true(sum(off) > 0)

    # A model built by this version says nothing.
    reset <- function() rm(list = ls(envir = ctsem:::.ct_legacy_covtransform),
      envir = ctsem:::.ct_legacy_covtransform)
    reset()
    expect_silent(ctsem:::.ctCheckLegacyCovTransformModel(m$pars))

    # One carrying the old expression is reported, once.
    legacy <- m
    legacy$pars$transform[off] <- ctsem:::.CT_LEGACY_COR_TRANSFORM
    reset()
    expect_warning(ctsem:::.ctCheckLegacyCovTransformModel(legacy$pars),
      'applied twice')
    expect_silent(ctsem:::.ctCheckLegacyCovTransformModel(legacy$pars))

    # A hand-written transform on the same cell is not second-guessed.
    custom <- m
    custom$pars$transform[off] <- 'tanh(param)'
    reset()
    expect_silent(ctsem:::.ctCheckLegacyCovTransformModel(custom$pars))

    # And the stored-standata shape, which is where a saved stan fit carries
    # it: transform code 3 on an off-diagonal of matrix 4, 5 or 8. Columns are
    # positional -- 1 row, 2 col, 4 transform, 7 matrix.
    p <- suppressWarnings(suppressMessages(ctFit(datalong = covmat_data(),
      model = m, fit = FALSE, cores = 1L, verbose = 0L)))
    sdat <- p$standata
    ms <- sdat$matsetup
    soff <- ms[, 7] %in% c(4, 5, 8) & ms[, 1] != ms[, 2] & ms[, 3] > 0
    expect_true(sum(soff) > 0)
    reset()
    expect_silent(ctsem:::.ctCheckLegacyCovTransform(sdat))
    sdat$matsetup[soff, 4] <- 3L
    reset()
    expect_warning(ctsem:::.ctCheckLegacyCovTransform(sdat), 'applied twice')
    reset()
  })
}
