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

  # The message is the only thing that tells a user their model lost population
  # parameters, so what it claims has to be right. An effect regressed because
  # its spread is unidentified and one regressed to meet a requested rank are
  # different claims: the first loses nothing, the second is an approximation.
  # A first version explained every regressed effect with the identification
  # reason and so asserted that a DRIFT effect "varies only in DIFFUSION /
  # MANIFESTVAR", which is false.
  test_that('the poprank message separates unidentified from approximated', {
    pars <- prepared_pars(poprank_model6())

    auto <- ctsem:::.ctPopRegressionSpec(pars, 'auto')
    auto$coefficients <- data.frame(param = auto$regressed,
      stringsAsFactors = FALSE)
    message_auto <- ctsem:::.ctPopRegressionMessage(auto)
    expect_match(message_auto, 'rank 3 of 6', fixed = TRUE)
    expect_match(message_auto, 'df1, df2, df3 vary', fixed = TRUE)
    expect_match(message_auto, 'cannot see their own spread', fixed = TRUE)
    # Nothing was approximated, so there is no approximation clause and no
    # effect is reported as demoted. Not asserted by searching for 'dr1': the
    # DRIFT effects are the *basis* here, so they are named legitimately in
    # "a regression on dr1, dr2, dr3".
    expect_false(grepl('approximation', message_auto, fixed = TRUE))
    expect_false(grepl('also regressed', message_auto, fixed = TRUE))
    expect_match(message_auto, 'regression on dr1, dr2, dr3', fixed = TRUE)
    expect_match(message_auto, 'poprank=NA', fixed = TRUE)

    one <- ctsem:::.ctPopRegressionSpec(pars, 1L)
    one$coefficients <- data.frame(param = one$regressed,
      stringsAsFactors = FALSE)
    message_one <- ctsem:::.ctPopRegressionMessage(one)
    expect_match(message_one, 'rank 1 of 6', fixed = TRUE)
    # the variance-cell effects keep the identification reason ...
    expect_match(message_one, 'df1, df2, df3 vary', fixed = TRUE)
    # ... and the demoted DRIFT effects are reported as the approximation
    expect_match(message_one, 'dr2, dr3 are also regressed', fixed = TRUE)
    expect_match(message_one, 'below the 3 this model identifies', fixed = TRUE)
    expect_match(message_one, 'approximation', fixed = TRUE)
  })

  # Under laplace and 'none' there are no carrier states, so the basis effects
  # move into PARS and the regressed cells reference them as parameters. No
  # engine change: the laplace spec builds re_index from the parameter table's
  # indvarying entries, so a PARS row is eligible exactly as the original cell
  # was, and the route simply sees fewer random effects.
  test_that('poprank restricts the laplace and none routes too, when asked', {
    dat <- poprank_data()
    prepared <- function(iop, ...) {
      suppressWarnings(suppressMessages(ctFit(datalong = dat,
        model = poprank_model(), backend = 'julia', fit = FALSE,
        intoverpop = iop, cores = 1L, ...)))
    }
    for (iop in c('laplace', 'augmented')) {
      full <- prepared(iop, poprank = NA)
      auto <- prepared(iop, poprank = 'auto')
      expect_equal(ctsem:::.ctBackendNpar(full), 6L)
      expect_equal(ctsem:::.ctBackendNpar(auto), 5L)
      names <- ctsem:::.ctBackendRawParameterNames(list(model_spec = auto),
        ctsem:::.ctBackendNpar(auto))
      expect_true('beta_df11_dr11' %in% names)
      expect_false(any(c('popsd_df11', 'rawcor_df11__dr11') %in% names))
    }
  })

  # And it is NOT the default off the augmented route. On the augmented route
  # 'auto' removes coordinates that cannot be identified, so it costs nothing.
  # Under laplace those coordinates are identified: measured on a 250 x 50
  # design the same restriction cost 48 log likelihood units and took the fit to
  # the boundary, basis sd to zero with the coefficient at -612. A default doing
  # that to someone who chose laplace *because* it identifies these things would
  # be indefensible.
  test_that('the poprank default does not restrict a laplace fit', {
    dat <- poprank_data()
    auto <- suppressWarnings(suppressMessages(ctFit(datalong = dat,
      model = poprank_model(), backend = 'julia', fit = FALSE,
      intoverpop = 'laplace', cores = 1L)))
    expect_null(auto$model$popregression)
    expect_equal(ctsem:::.ctBackendNpar(auto), 6L)
    # the augmented default does restrict, which is the contrast
    augmented <- suppressWarnings(suppressMessages(ctFit(datalong = dat,
      model = poprank_model(), backend = 'julia', fit = FALSE,
      intoverpop = 'augmented', cores = 1L)))
    expect_equal(ctsem:::.ctBackendNpar(augmented), 5L)
  })

  # The message must not carry the augmented route's justification onto a route
  # where it is false. A first version told a laplace user that df11 "varies
  # only in DIFFUSION / MANIFESTVAR, where the augmented filter cannot see its
  # own spread" -- on the one route that does identify it.
  test_that('the poprank message says which claim it is making', {
    pars <- prepared_pars(poprank_model())
    spec <- ctsem:::.ctPopRegressionSpec(pars, 'auto')
    spec$coefficients <- data.frame(param = spec$regressed,
      stringsAsFactors = FALSE)

    spec$route <- 'augmented'
    augmented <- ctsem:::.ctPopRegressionMessage(spec)
    expect_match(augmented, 'cannot see its own spread', fixed = TRUE)
    expect_match(augmented, "intoverpop='laplace' identifies it", fixed = TRUE)

    spec$route <- 'parameters'
    parameters <- ctsem:::.ctPopRegressionMessage(spec)
    expect_false(grepl('cannot see', parameters, fixed = TRUE))
    expect_match(parameters, 'those coordinates are identified', fixed = TRUE)
    expect_match(parameters, 'approximation for parsimony or speed', fixed = TRUE)
  })

  # The rank may be stated on the model, which is the idiom indvarying already
  # uses. An argument at the call site wins over it.
  test_that('poprank can be stated on the model, and the argument wins', {
    dat <- poprank_data()
    m <- poprank_model()
    m$poprank <- 'auto'
    onmodel <- suppressWarnings(suppressMessages(ctFit(datalong = dat, model = m,
      backend = 'julia', fit = FALSE, intoverpop = 'augmented', cores = 1L)))
    expect_equal(ctsem:::.ctBackendNpar(onmodel), 5L)

    overridden <- suppressWarnings(suppressMessages(ctFit(datalong = dat,
      model = m, backend = 'julia', fit = FALSE, intoverpop = 'augmented',
      poprank = NA, cores = 1L)))
    expect_equal(ctsem:::.ctBackendNpar(overridden), 6L)
  })

  # What POPCOV can and cannot say once the rank is reduced. The freely
  # parameterised part keeps today's sd and correlation coordinates, so
  # statements there work untouched; a zero *variance* has one clear meaning
  # under any rank and is honoured by leaving that effect out of the split; and
  # a non-zero statement about a regressed effect, or a zero covariance -- which
  # is a linear constraint across a whole row of coefficients rather than a cell
  # -- is refused rather than dropped.
  test_that('POPCOV statements work where they can be honoured and are refused where not', {
    m <- poprank_model6()
    pars <- prepared_pars(m)
    m$pars <- pars
    fresh <- function() { out <- m; out[['POPCOV']] <- ctsem:::.ctModelPopCov(pars); out }
    split <- function(model) ctsem:::.ctPopRegressionSpec(prepared_pars(model),
      'auto', explicit = TRUE, model = model)

    base <- split(fresh())
    expect_equal(base$basis, c('dr1', 'dr2', 'dr3'))
    expect_equal(base$regressed, c('df1', 'df2', 'df3'))

    # a zero correlation between two basis effects: the freely parameterised part
    z <- fresh(); z[['POPCOV']]['dr2', 'dr1'] <- 0
    expect_equal(split(z)$regressed, c('df1', 'df2', 'df3'))

    # a fixed sd on a basis effect
    f <- fresh(); f[['POPCOV']]['dr1', 'dr1'] <- 0.3
    expect_equal(split(f)$basis, c('dr1', 'dr2', 'dr3'))

    # and any statement about a regressed effect is refused, zero or not.
    # A zero on the diagonal is not special-cased: fixing a population sd to
    # zero is not a sensible thing to state -- `indvarying = FALSE` is how a
    # parameter is made non-varying -- so it gets no path of its own.
    nz <- fresh(); nz[['POPCOV']]['df1', 'df1'] <- 0.3
    expect_error(split(nz), 'poprank would drop what POPCOV states')
    zd <- fresh(); zd[['POPCOV']]['df1', 'df1'] <- 0
    expect_error(split(zd), 'poprank would drop what POPCOV states')
    zc <- fresh(); zc[['POPCOV']]['df1', 'dr1'] <- 0
    expect_error(split(zc), 'poprank would drop what POPCOV states')
  })

  # What a fixed POPCOV entry means. Pinned because the documentation in
  # R/ctModelPopCov.R now states these numbers, and a comment is not evidence:
  # if `constraincorsqrt1()` or the placement changes, this is what notices.
  #
  # The diagonal is the population sd it says it is. The off-diagonal is the
  # coordinate that map consumes, so the correlation that comes out is further
  # from zero than the coordinate written in -- except at zero, which is exact
  # and is the case the surface's own example relies on.
  test_that('a fixed POPCOV entry means what R/ctModelPopCov.R says it means', {
    set.seed(4); nsub <- 120L; nt <- 6L
    o <- vector('list', nsub)
    for (i in seq_len(nsub)) {
      eta <- numeric(nt); eta[1] <- rnorm(1)
      for (t in 2:nt) eta[t] <- exp(-.5) * eta[t - 1] + rnorm(1, 0, .4)
      o[[i]] <- data.frame(id = i, time = seq_len(nt) - 1,
        Y1 = eta + rnorm(1, 0, .5) + rnorm(nt, 0, .3),
        Y2 = eta + rnorm(1, 0, .5) + rnorm(nt, 0, .3))
    }
    dat <- do.call(rbind, o)
    mk <- function() {
      m <- suppressMessages(ctModel(type = 'ct', n.latent = 1, n.manifest = 2,
        LAMBDA = matrix(c(1, 1), 2, 1), manifestNames = c('Y1', 'Y2'),
        latentNames = 'eta', DRIFT = matrix(-0.5), DIFFUSION = matrix(0.4),
        CINT = matrix(0), T0MEANS = matrix(0), T0VAR = matrix(1),
        MANIFESTMEANS = matrix(c('mm1|param', 'mm2|param'), 2, 1),
        MANIFESTVAR = matrix(c('mv1|log1p_exp(param)', 0, 0,
          'mv2|log1p_exp(param)'), 2, 2)))
      m$pars$indvarying <- FALSE
      m$pars$indvarying[m$pars$matrix %in% 'MANIFESTMEANS'] <- TRUE
      m
    }
    withcell <- function(row, col, value) {
      m <- mk()
      p <- ctsem:::ctModelStatesAndPARS(ctsem:::ctModel0DRIFT(m, TRUE)$pars,
        statenames = 'eta', tdprednames = NULL)
      m$pars <- p
      m[['POPCOV']] <- ctsem:::.ctModelPopCov(p)
      nms <- rownames(m[['POPCOV']])
      m[['POPCOV']][nms[row], nms[col]] <- value
      set.seed(77)
      f <- suppressWarnings(suppressMessages(ctFit(datalong = dat, model = m,
        backend = 'julia', intoverpop = 'augmented', poprank = NA, cores = 1L,
        verbose = 0L, optimcontrol = list(estonly = TRUE))))
      ctsem:::.ctBackendRawPopCov(f)[[1L]]$cov
    }
    correlation <- function(cov) cov[2, 1] / sqrt(cov[1, 1] * cov[2, 2])

    # the diagonal is the sd, exactly, for these linear transforms
    expect_equal(sqrt(withcell(1, 1, 0.3)[1, 1]), 0.3, tolerance = 1e-4)
    expect_equal(sqrt(withcell(1, 1, 0.8)[1, 1]), 0.8, tolerance = 1e-4)

    # zero off the diagonal is uncorrelated, exactly
    expect_equal(correlation(withcell(2, 1, 0)), 0, tolerance = 1e-6)

    # and a non-zero coordinate gives a correlation further from zero than
    # itself -- the numbers the documentation quotes
    expect_equal(correlation(withcell(2, 1, 0.3)), 0.539, tolerance = 2e-3)
    expect_equal(correlation(withcell(2, 1, 0.5)), 0.788, tolerance = 2e-3)
    # Monotone and sign-preserving, which is why it reads like a correlation --
    # but not symmetric in sign, which is another way it is not one. The same
    # coordinate magnitude gives a different correlation magnitude either side
    # of zero, because constraincorsqrt1()'s row scale carries an |s| - s term.
    expect_equal(correlation(withcell(2, 1, -0.5)), -0.769, tolerance = 2e-3)
    expect_lt(abs(correlation(withcell(2, 1, -0.5))),
      abs(correlation(withcell(2, 1, 0.5))))
  })

  # ctIdentify deliberately assesses the *unrestricted* covariance, whatever
  # ctFit's poprank default is, and this pins that rather than leaving it to be
  # "fixed" later. Its job is to say what the data identifies -- nweak on the
  # full-rank model is exactly why poprank exists, and reporting the reduced
  # model's 0 instead would hide the finding. Routing it through
  # ctFit(fit = FALSE) to pick up poprank was tried and reverted: it made
  # ctIdentify consume random numbers (it is documented as leaving the RNG
  # alone), broke intoverpop='laplace' through it, and changed an error
  # message, for a diagnostic that was already answering the right question.
  test_that('ctIdentify assesses the unrestricted covariance, not the reduced one', {
    dat <- poprank_data(nsub = 60L)
    id <- suppressWarnings(suppressMessages(ctIdentify(dat, poprank_model(),
      nstart = 1L, cores = 1L, verbose = 0L)))
    # three population parameters, and the flat direction poprank removes
    expect_equal(id$npar, 6L)
    expect_gte(id$nweak, 1L)
    expect_false('poprank' %in% names(formals(ctIdentify)))
  })

  # POPCOV is the specification surface for the population covariance, and the
  # augmentation reads it per varying parameter: a number fixes an sd or a
  # correlation. A regressed effect has neither of its own, so anything stated
  # about one would be silently dropped -- the single outcome this feature
  # exists to prevent. Asked for, refused by name; defaulted, the user's own
  # specification wins and the rank is simply not applied.
  test_that('poprank refuses to drop a POPCOV statement about a regressed effect', {
    m <- poprank_model()
    pars <- prepared_pars(m)
    m$pars <- pars
    # a fresh POPCOV over the two varying parameters, then a fixed sd for the
    # one poprank would regress
    m[['POPCOV']] <- ctsem:::.ctModelPopCov(pars)
    expect_true('df11' %in% rownames(m[['POPCOV']]))
    m[['POPCOV']]['df11', 'df11'] <- 0.3

    conflicts <- ctsem:::.ctPopRegressionPopCovConflicts(m, 'df11')
    expect_length(conflicts, 1L)
    expect_match(conflicts, "POPCOV['df11', 'df11'] = 0.3", fixed = TRUE)

    expect_error(ctsem:::.ctPopRegressionSpec(pars, 'auto', explicit = TRUE,
      model = m), 'poprank would drop what POPCOV states')
    expect_null(ctsem:::.ctPopRegressionSpec(pars, 'auto', explicit = FALSE,
      model = m))

    # a statement about the basis effect is fine -- it keeps its own spread
    m2 <- m
    m2[['POPCOV']] <- ctsem:::.ctModelPopCov(pars)
    m2[['POPCOV']]['dr11', 'dr11'] <- 0.3
    expect_length(ctsem:::.ctPopRegressionPopCovConflicts(m2, 'df11'), 0L)
    spec <- ctsem:::.ctPopRegressionSpec(pars, 'auto', explicit = TRUE, model = m2)
    expect_equal(spec$regressed, 'df11')

    # and an untouched POPCOV states nothing, so it cannot conflict
    m3 <- m
    m3[['POPCOV']] <- ctsem:::.ctModelPopCov(pars)
    expect_length(ctsem:::.ctPopRegressionPopCovConflicts(m3, 'df11'), 0L)
  })

  test_that('poprank refuses a model in which nothing reaches the observation mean', {
    m <- poprank_model()
    m$pars$indvarying <- FALSE
    m$pars$indvarying[m$pars$matrix %in% 'DIFFUSION'] <- TRUE
    expect_error(ctsem:::.ctPopRegressionSpec(prepared_pars(m), 'auto'),
      'reaches the observation mean')
  })

  # Explicitly asking for a rank where it cannot apply is an error; the default
  # is simply not applied, which is what lets 'auto' be the default at all.
  # Laplace is no longer among the refusals -- it is supported, just not by
  # default -- so only the backend is.
  test_that('an explicitly requested poprank is refused on stan', {
    dat <- poprank_data()
    expect_error(suppressWarnings(ctFit(datalong = dat, model = poprank_model(),
      backend = 'stan', fit = FALSE, intoverpop = 'augmented', poprank = 'auto')),
      "requires backend='julia'")
  })

  test_that('the default poprank is silently inapplicable where it cannot be used', {
    dat <- poprank_data()
    expect_type(suppressWarnings(suppressMessages(ctFit(datalong = dat,
      model = poprank_model(), backend = 'stan', fit = FALSE,
      intoverpop = 'augmented'))), 'list')
    laplace <- suppressWarnings(suppressMessages(ctFit(datalong = dat,
      model = poprank_model(), backend = 'julia', fit = FALSE,
      intoverpop = 'laplace', cores = 1L)))
    expect_null(laplace$model$popregression)
  })

  # A model where nothing reaches the observation mean has no basis to regress
  # on. Asked for, that is an error; defaulted, the model is left alone rather
  # than a call that used to run becoming a failure.
  test_that('the default poprank leaves a model with no mean-affecting effect alone', {
    dat <- poprank_data()
    m <- poprank_model()
    m$pars$indvarying <- FALSE
    m$pars$indvarying[m$pars$matrix %in% 'DIFFUSION'] <- TRUE
    expect_null(ctsem:::.ctPopRegressionSpec(prepared_pars(m), 'auto',
      explicit = FALSE))
    expect_error(ctsem:::.ctPopRegressionSpec(prepared_pars(m), 'auto',
      explicit = TRUE), 'reaches the observation mean')
    spec <- suppressWarnings(suppressMessages(ctFit(datalong = dat, model = m,
      backend = 'julia', fit = FALSE, intoverpop = 'augmented', cores = 1L)))
    expect_null(spec$model$popregression)
  })

  # 'auto' is the default, so the case where it changes nothing has to be
  # verified rather than assumed: with both effects mean-affecting the rank is
  # already full and the model must come out exactly as poprank=NA does.
  test_that('poprank auto is a no-op when every varying parameter reaches the mean', {
    dat <- poprank_data()
    m <- poprank_model()
    m$pars$indvarying <- FALSE
    m$pars$indvarying[m$pars$matrix %in% c('DRIFT', 'CINT', 'T0MEANS') &
        !is.na(m$pars$param)] <- TRUE
    # only DRIFT is a free parameter here, so k = 1 and nothing is regressed
    spec <- ctsem:::.ctPopRegressionSpec(prepared_pars(m), 'auto')
    expect_equal(spec$rank, 1L)
    expect_length(spec$regressed, 0L)
    auto <- suppressWarnings(suppressMessages(ctFit(datalong = dat, model = m,
      backend = 'julia', fit = FALSE, intoverpop = 'augmented', cores = 1L)))
    none <- suppressWarnings(suppressMessages(ctFit(datalong = dat, model = m,
      backend = 'julia', fit = FALSE, intoverpop = 'augmented',
      poprank = NA, cores = 1L)))
    expect_equal(ctsem:::.ctBackendNpar(auto), ctsem:::.ctBackendNpar(none))
    expect_null(auto$model$popregression)
  })

  # The fit-free guard: cheap, and it is what catches a rewrite that dropped a
  # parameter or left one behind. Full rank has 6 (mv, dr11, df11, two
  # population sds, one correlation); the regression form replaces the second
  # sd and the correlation with a single coefficient, and needs one carrier
  # state rather than two.
  test_that('poprank removes the unidentified coordinates and one carrier state', {
    dat <- poprank_data()
    # poprank defaults to 'auto', so the unrestricted arm has to ask for NA.
    full <- suppressWarnings(suppressMessages(ctFit(datalong = dat,
      model = poprank_model(), backend = 'julia', fit = FALSE,
      intoverpop = 'augmented', poprank = NA, cores = 1L)))
    auto <- suppressWarnings(suppressMessages(ctFit(datalong = dat,
      model = poprank_model(), backend = 'julia', fit = FALSE,
      intoverpop = 'augmented', poprank = 'auto', cores = 1L)))

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

  # What a reduced-rank fit shows a user. The point of the feature is that the
  # interesting quantities are the same ones as always -- how much each
  # parameter varies and how the spreads go together -- so the tables keep their
  # shape and the *structure* is carried by a note. The regression coefficients
  # are the mechanism and stay out of the way: not a printed section, and not a
  # row in the population means beside the model's own parameters.
  # 100 subjects rather than this file's 40, deliberately. At 40 the rank-1
  # optimum on this data sits at the boundary -- the basis spread goes to zero
  # with the coefficient going to infinity, holding their product -- so
  # `popsd` is legitimately 0 there and an assertion that the spreads are
  # positive would be asserting something false. Verified that -384.4209 is the
  # genuine rank-1 optimum and not an optimiser failure: a grid over sensible
  # (basis sd, coefficient) values, then optimised, reaches only -384.467. Full
  # rank gets -383.846 on the same data, so at 40 subjects the restriction
  # costs 0.575 -- it refuses to spend an unidentified parameter on noise, which
  # is the point, but it is not the "costs nothing" of the well-determined case.
  # ctsem already reports that fit properly: non-convergence, and the
  # identifiability warning naming popsd_dr11 and beta_df11_dr11 as not
  # estimable with NA widths. No extra warning was added for it.
  test_that('a poprank fit reports the spreads, not the mechanism', {
    dat <- poprank_data(nsub = 100L)
    set.seed(303)
    f <- suppressWarnings(suppressMessages(ctFit(datalong = dat,
      model = poprank_model(), backend = 'julia', intoverpop = 'augmented',
      cores = 1L, verbose = 0L)))
    s <- suppressWarnings(suppressMessages(summary(f, parmatrices = FALSE,
      priorcheck = FALSE, residualcov = FALSE)))

    # every varying parameter has a spread, including the regressed one
    expect_setequal(rownames(s$popsd), c('dr11', 'df11'))
    expect_true(all(s$popsd[, 'mean'] > 0))
    # and a correlation between them
    expect_equal(rownames(s$rawpopcorr), 'df11__dr11')

    # the note says the dimension structure
    expect_false(is.null(s$popsdNote))
    expect_match(s$popsdNote, '1 dimension, not 2', fixed = TRUE)
    expect_match(s$popsdNote, 'no variation independent of dr11', fixed = TRUE)
    expect_match(s$popsdNote, 'DIFFUSION / MANIFESTVAR', fixed = TRUE)

    # the mechanism is not in the way
    expect_null(s$popregression)
    expect_false(any(grepl('^beta_', rownames(s$popmeans))))
    # ... while the parameter count still counts it
    expect_equal(length(f$estimate$raw), 5L)
    # ... and it is still reachable for anyone who wants it
    coefficients <- ctsem:::.ctBackendPopRegressionTable(f,
      ctsem:::.ctBackendSpec(f), ctsem:::.ctBackendRawSamples(f))
    expect_equal(coefficients$param, 'df11')
    expect_equal(coefficients$on, 'dr11')

    # the resolved rank is recorded as a number, not as the argument
    expect_equal(f$args$resolved$poprank, 1L)

    # and a figure says it too, since a figure outlives the fit message
    tex <- suppressWarnings(suppressMessages(ctModelLatex(f, compile = FALSE,
      open = FALSE, equationonly = TRUE)))
    expect_match(paste(tex, collapse = ' '),
      'Individual differences have 1 dimension of 2', fixed = TRUE)
  })

  # A regressed effect varies by subject without a carrier state of its own, so
  # the enumeration ctSubjectPars is built on could not see it and it was
  # silently absent -- the failure mode that looks identical to "not
  # applicable" from outside.
  #
  # Nothing here compares the two fits' per-subject *values*, and the reason is
  # worth recording rather than discovering again. Under full rank the position
  # on the ridge is arbitrary, and it perturbs the per-subject estimates as well
  # as the population sd -- not only for the variance-cell effect whose carrier
  # state the filter never updates, but for the identified effect too, because
  # the two states are correlated and the arbitrary one changes the gain.
  # Measured on this dataset, at likelihoods agreeing to 1e-10: per-subject
  # values differing by up to 0.05, with Spearman correlations between the two
  # fits of 0.93 for `dr11` and 0.89 for `df11`. So neither equality nor a
  # rank-correlation threshold is an invariant here; picking a threshold that
  # passes would be a test that discriminates nothing. What is asserted below is
  # what actually holds: the column exists, it varies by subject, and under
  # rank 1 the regressed effect is a monotone function of the basis by
  # construction.
  test_that('ctSubjectPars carries a regressed effect and it varies by subject', {
    dat <- poprank_data()
    fitboth <- function(poprank) {
      set.seed(303)
      args <- list(datalong = dat, model = poprank_model(), backend = 'julia',
        intoverpop = 'augmented', cores = 1L, verbose = 0L)
      args$poprank <- poprank
      suppressWarnings(suppressMessages(do.call(ctFit, args)))
    }
    auto <- fitboth('auto')
    full <- fitboth(NA)

    pauto <- suppressWarnings(suppressMessages(ctSubjectPars(auto)))
    pfull <- suppressWarnings(suppressMessages(ctSubjectPars(full)))

    expect_setequal(dimnames(pauto)$param, c('dr11', 'df11'))
    expect_setequal(dimnames(pfull)$param, c('dr11', 'df11'))
    expect_equal(dim(pauto), dim(pfull))

    for (p in c('dr11', 'df11')) {
      # genuinely per subject, not a constant column, on both routes
      expect_gt(stats::sd(apply(pauto[, , p, drop = FALSE], 2L, mean)), 0)
      expect_gt(stats::sd(apply(pfull[, , p, drop = FALSE], 2L, mean)), 0)
    }

    # And under rank 1 the regressed effect is a monotone function of the basis
    # effect by construction, which is the structure rather than an estimate.
    values <- apply(pauto, c(2, 3), mean)
    expect_equal(abs(stats::cor(values[, 'dr11'], values[, 'df11'],
      method = 'spearman')), 1)
  })

  # A TI effect shifts a subject's raw parameter value, which under the
  # augmented route is the population mean rather than the random deviation.
  # So a TI effect on a regressed effect follows its mean parameter, and must
  # not be left on the cell that becomes an expression -- a flag there makes
  # `.ctJuliaTIEffects()` mint a coefficient with nothing to attach it to.
  test_that('a TI-predictor effect on a regressed effect follows its mean', {
    m <- poprank_model()
    m$pars$TI1_effect <- FALSE
    m$TIpredNames <- 'TI1'
    m$n.TIpred <- 1L
    m$pars$TI1_effect[!is.na(m$pars$param) & m$pars$param %in% 'df11'] <- TRUE
    pars <- prepared_pars(m)
    m$pars <- pars
    spec <- ctsem:::.ctPopRegressionSpec(pars, 'auto')
    m <- ctsem:::.ctPopRegressionDemote(m, spec)
    m <- ctsem:::.ctModelIntOverPop(m)
    m <- ctsem:::.ctPopRegressionRewrite(m, spec)

    mean <- m$pars$matrix %in% 'PARS' & m$pars$param %in% 'df11'
    expect_true(any(mean))
    expect_true(all(m$pars$TI1_effect[mean]))
    # and nowhere on an expression cell
    expression <- !is.na(m$pars$param) & grepl('[', m$pars$param, fixed = TRUE)
    expect_false(any(m$pars$TI1_effect[expression]))
    # the coefficients themselves take no TI effect
    beta <- m$pars$param %in% 'beta_df11_dr11'
    expect_true(any(beta))
    expect_false(any(m$pars$TI1_effect[beta]))
  })

}
