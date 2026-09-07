if (identical(Sys.getenv('NOT_CRAN'), 'true')) {

  library(ctsem)
  library(testthat)

  # A one-latent model with a DRIFT random effect (mean-affecting, so
  # identified) and a DIFFUSION random effect (a variance cell, so only its
  # covariance with the DRIFT effect is identified).
  poprank_model <- function() {
    m <- suppressMessages(ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
      LAMBDA = matrix(1), manifestNames = 'Y1', latentNames = 'eta',
      DRIFT = matrix('dr11|-log1p_exp(param)'),
      DIFFUSION = matrix('df11|log1p_exp(param)'),
      CINT = matrix(0), MANIFESTMEANS = matrix(0),
      MANIFESTVAR = matrix('mv|log1p_exp(param)'),
      T0MEANS = matrix(0), T0VAR = matrix(1)))
    m$pars$indvarying <- FALSE
    m$pars$indvarying[m$pars$matrix %in% c('DRIFT', 'DIFFUSION')] <- TRUE
    m
  }

  # Three latents, three DRIFT and three DIFFUSION random effects: m = 3, v = 3.
  poprank_model6 <- function() {
    nl <- 3L
    DR <- matrix(0, nl, nl); diag(DR) <- paste0('dr', 1:nl, '|-log1p_exp(param)')
    DF <- matrix(0, nl, nl); diag(DF) <- paste0('df', 1:nl, '|log1p_exp(param)')
    MV <- matrix(0, nl, nl); diag(MV) <- paste0('mv', 1:nl, '|log1p_exp(param)')
    m <- suppressMessages(ctModel(type = 'ct', n.latent = nl, n.manifest = nl,
      LAMBDA = diag(nl), manifestNames = paste0('Y', 1:nl),
      latentNames = paste0('eta', 1:nl), DRIFT = DR, DIFFUSION = DF,
      MANIFESTVAR = MV, CINT = matrix(0, nl, 1), MANIFESTMEANS = matrix(0, nl, 1),
      T0MEANS = matrix(0, nl, 1), T0VAR = diag(1, nl)))
    m$pars$indvarying <- FALSE
    m$pars$indvarying[m$pars$matrix %in% c('DRIFT', 'DIFFUSION') &
        !is.na(m$pars$param)] <- TRUE
    m
  }

  prepared_pars <- function(m) {
    m <- ctsem:::ctModel0DRIFT(m, m$continuoustime)
    ctsem:::ctModelStatesAndPARS(m$pars, statenames = m$latentNames,
      tdprednames = m$TDpredNames)
  }

  poprank_data <- function(nsub = 40L, nt = 8L, seed = 21L) {
    set.seed(seed)
    z <- matrix(rnorm(nsub * 2), nsub, 2) %*% chol(matrix(c(1, .7, .7, 1), 2, 2))
    out <- vector('list', nsub)
    for (i in seq_len(nsub)) {
      a <- -log1p(exp(z[i, 1] * .6)); q <- log1p(exp(z[i, 2] * .4))
      e <- exp(a); eta <- numeric(nt); eta[1] <- rnorm(1)
      for (t in 2:nt) eta[t] <- e * eta[t - 1] +
        rnorm(1, 0, sqrt(q^2 / (-2 * a) * (1 - e^2)))
      out[[i]] <- data.frame(id = i, time = seq_len(nt) - 1,
        Y1 = eta + rnorm(nt, 0, sqrt(.2)))
    }
    do.call(rbind, out)
  }

  test_that('a DRIFT random effect reaches the observation mean and a DIFFUSION one does not', {
    roles <- ctsem:::.ctPopEffectRoles(prepared_pars(poprank_model()))
    expect_equal(sort(roles$param), c('df11', 'dr11'))
    expect_true(roles$mean[roles$param == 'dr11'])
    expect_false(roles$mean[roles$param == 'df11'])
  })

  # This is the case a pattern match over `$pars$matrix` gets wrong, and the
  # reason the classifier follows PARS references rather than reading the
  # matrix column: `df11` drives DIFFUSION *and* appears inside the DRIFT
  # expression, so it does reach the observation mean and is identified, while
  # having no DRIFT row of its own.
  test_that('a variance-cell parameter used in a mean expression counts as mean-affecting', {
    m <- suppressMessages(ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
      LAMBDA = matrix(1), manifestNames = 'Y1', latentNames = 'eta',
      PARS = matrix(c('df11', 'dr11'), 2, 1),
      DRIFT = matrix('-log1p_exp(dr11 + 0.1 * df11)'),
      DIFFUSION = matrix('log1p_exp(df11)'),
      CINT = matrix(0), MANIFESTMEANS = matrix(0),
      MANIFESTVAR = matrix('mv|log1p_exp(param)'),
      T0MEANS = matrix(0), T0VAR = matrix(1)))
    m$pars$indvarying <- FALSE
    m$pars$indvarying[!is.na(m$pars$param) &
        m$pars$param %in% c('df11', 'dr11')] <- TRUE
    expect_true(all(ctsem:::.ctPopEffectRoles(prepared_pars(m))$mean))

    # And the same shape without the DRIFT reference is variance-only again,
    # so the test above is not passing for want of discrimination.
    m$pars$param[m$pars$matrix %in% 'DRIFT'] <- '-log1p_exp(dr11)'
    roles <- ctsem:::.ctPopEffectRoles(prepared_pars(m))
    expect_false(roles$mean[roles$param == 'df11'])
    expect_true(roles$mean[roles$param == 'dr11'])
  })

  test_that('poprank resolves to the right basis and parameter count', {
    pars <- prepared_pars(poprank_model6())
    expect_null(ctsem:::.ctPopRegressionSpec(pars, NA))

    auto <- ctsem:::.ctPopRegressionSpec(pars, 'auto')
    expect_equal(auto$rank, 3L)
    expect_equal(auto$basis, c('dr1', 'dr2', 'dr3'))
    expect_equal(auto$regressed, c('df1', 'df2', 'df3'))
    # m(m+1)/2 + v*m against the full-rank k(k+1)/2 = 21
    expect_equal(auto$npar, 15L)
    expect_false(auto$approximate)

    expect_equal(ctsem:::.ctPopRegressionSpec(pars, 2L)$npar, 11L)
    expect_equal(ctsem:::.ctPopRegressionSpec(pars, 1L)$npar, 6L)
    expect_true(ctsem:::.ctPopRegressionSpec(pars, 2L)$approximate)
    # An explicit rank at the identified dimension is not an approximation.
    expect_false(ctsem:::.ctPopRegressionSpec(pars, 3L)$approximate)
    expect_error(ctsem:::.ctPopRegressionSpec(pars, 7L), 'only 6')
  })

  test_that('poprank refuses a model in which nothing reaches the observation mean', {
    m <- poprank_model()
    m$pars$indvarying <- FALSE
    m$pars$indvarying[m$pars$matrix %in% 'DIFFUSION'] <- TRUE
    expect_error(ctsem:::.ctPopRegressionSpec(prepared_pars(m), 'auto'),
      'reaches the observation mean')
  })

  test_that('poprank is refused on stan and under intoverpop=laplace', {
    dat <- poprank_data()
    expect_error(suppressWarnings(ctFit(datalong = dat, model = poprank_model(),
      backend = 'stan', fit = FALSE, intoverpop = 'augmented', poprank = 'auto')),
      "requires backend='julia'")
    expect_error(suppressWarnings(ctFit(datalong = dat, model = poprank_model(),
      backend = 'julia', fit = FALSE, intoverpop = 'laplace', poprank = 'auto')),
      'augmented')
  })

  # The fit-free guard: cheap, and it is what catches a rewrite that dropped a
  # parameter or left one behind. Full rank has 6 (mv, dr11, df11, two
  # population sds, one correlation); the regression form replaces the second
  # sd and the correlation with a single coefficient, and needs one carrier
  # state rather than two.
  test_that('poprank removes the unidentified coordinates and one carrier state', {
    dat <- poprank_data()
    full <- suppressWarnings(ctFit(datalong = dat, model = poprank_model(),
      backend = 'julia', fit = FALSE, intoverpop = 'augmented', cores = 1L))
    auto <- suppressWarnings(ctFit(datalong = dat, model = poprank_model(),
      backend = 'julia', fit = FALSE, intoverpop = 'augmented',
      poprank = 'auto', cores = 1L))

    nfull <- ctsem:::.ctBackendNpar(full)
    nauto <- ctsem:::.ctBackendNpar(auto)
    expect_equal(nfull, 6L)
    expect_equal(nauto, 5L)

    namesfull <- ctsem:::.ctBackendRawParameterNames(list(model_spec = full), nfull)
    namesauto <- ctsem:::.ctBackendRawParameterNames(list(model_spec = auto), nauto)
    expect_true(all(c('popsd_dr11', 'popsd_df11', 'rawcor_df11__dr11') %in% namesfull))
    expect_true('beta_df11_dr11' %in% namesauto)
    # The coordinates the profile likelihood cannot distinguish are gone.
    expect_false(any(c('popsd_df11', 'rawcor_df11__dr11') %in% namesauto))
    # ... and the population mean of the regressed effect is still estimated.
    expect_true(all(c('dr11', 'df11', 'popsd_dr11') %in% namesauto))

    augdim <- function(spec) {
      table <- as.data.frame(spec$parameter_table)
      max(table$row[table$matrix == 'T0MEANS'])
    }
    expect_equal(augdim(full), 3L)
    expect_equal(augdim(auto), 2L)
  })

  test_that('poprank cannot yet be combined with TI-predictor effects on a regressed effect', {
    m <- poprank_model()
    m$pars$TI1_effect <- FALSE
    m$TIpredNames <- 'TI1'
    m$n.TIpred <- 1L
    m$pars$TI1_effect[!is.na(m$pars$param) & m$pars$param %in% 'df11'] <- TRUE
    spec <- ctsem:::.ctPopRegressionSpec(prepared_pars(m), 'auto')
    expect_error(ctsem:::.ctPopRegressionDemote(m, spec), 'TI-predictor')
  })

}
