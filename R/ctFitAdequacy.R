# Adequacy of a fitted model, as tests.
#
# Two kinds of family. Most are built on the one-step-ahead prediction errors:
# a row's errors, each variable's taken given the row's others and
# standardised (variable-labelled conditional residuals, from the SEM misfit
# work in ../papers/2025 SEM misfit), are standard normal and independent of
# every other row of the unit under the model, whatever the spacing -- the
# prediction-error decomposition is the time-series form of whitening. So any
# per-unit sum of products of them has mean zero, and is a conditional moment
# test (Newey 1985; Tauchen 1985): its estimation effect is removed by
# projecting each unit's contribution on the per-unit scores, and its variance
# is the spread of the projected contributions over units, as the score tests
# do (R/ctScoreTest.R). The rest compare statistics of the data with the same
# statistics of datasets generated from the fit, which uses the data twice
# and is conservative (Robins, van der Vaart & Ventura 2000); the lagged
# correlations of the data never rejected in 300 simulated datasets under any
# departure, and are no longer run unless asked for.
#
# A statistic's variance is the model's own -- the errors are independent over
# rows, so a sum of products over row pairs has the sum of their squares as
# its variance -- less what estimation took. Within a family, many correlated
# statistics: their largest |z| and their sum of squares, calibrated with the
# statistics' correlation, and combined with any exact p-values the family has by
# the Cauchy combination (Liu & Xie 2020), which holds under any dependence.
# The sum has the power against many small departures -- a flat residual
# autocorrelation -- and the max against a few large ones. Holm across
# families. Design and simulations in
# review/MISFIT-multivariate-literature-2026-10-10.md (CT-SEM).

#' Adequacy of a fitted model, as tests
#'
#' Tests a fitted model in families, each aimed at one way a continuous time
#' model can fall short, and gives each family a p-value. Most families use the
#' one-step-ahead prediction errors, which under the model are independent and
#' standard normal whatever the spacing of the observations; each variable's
#' error is taken given the others observed at the same time, so a result names
#' a variable.
#'
#' \describe{
#'   \item{residuals}{Correlations of the errors within units, for every pair
#'     of rows, by time gap (\code{gaps}): the residual autocorrelation and
#'     cross-correlation in continuous time. Slow structure the model lacks
#'     shows as correlation over the gaps it spans.}
#'   \item{latent}{The same for each latent process: the update each row makes
#'     to the predicted state, standardised, correlated with the same subject's
#'     updates at other times. One per process rather than one per item, so a
#'     misfit in the dynamics is not diluted by item-specific structure.}
#'   \item{volatility}{The same for the squared errors, per variable: variance
#'     that clusters in time.}
#'   \item{pairs}{Each pair of variables' errors in the same row against the
#'     correlation the model implies: a wrong measurement error covariance or
#'     loading.}
#'   \item{distribution}{Skewness and kurtosis of each variable's errors, and
#'     each row's squared whitened distance against its chi-square mean and
#'     variance, together with Anderson-Darling tests against the normal and
#'     chi-square. On a Laplace fit whose random effects reach the dynamics or
#'     the variances, each row's predictive covariance is a Gaussian
#'     approximation, and the row distances are left out: they would test the
#'     approximation, and in simulation rejected correct models.}
#'   \item{fit}{How unevenly units are fitted: the spread of each unit's
#'     summed distance against the chi-square's, and the worst unit.}
#'   \item{dependence}{The distance correlation of each variable's error, and
#'     of its square, with what came before -- the predicted latent states, the
#'     previous errors, the time interval -- by permutation: a conditional mean
#'     or variance the model gets wrong, in any form.}
#'   \item{effects}{On a single-level Laplace fit, whether the random effects'
#'     normal distribution is adequate: each subject's posterior probability,
#'     smoothed, that its standardised effect lies below a cut, at cuts along
#'     each axis, against the normal's. Under the model the posterior averages
#'     to the prior, so the differences have mean zero; a skewed or bimodal
#'     distribution, or a subgroup, moves them. These integrate the gradient
#'     function (Verbeke & Molenberghs, 2013; Drikvandi, Verbeke & Molenberghs,
#'     2017), which is too heavy tailed to test directly.}
#'   \item{scores}{Not run unless asked for. White's (1987) dynamic information
#'     matrix test: each data row's score, the derivative of its one-step
#'     log likelihood, for the dynamic and measurement parameters
#'     (\code{parameters}), correlated with the scores of the same subject's
#'     other rows, by time gap. Under the model a subject's row scores are
#'     uncorrelated whatever the departure, so this is a check for dynamic
#'     misspecification of any kind; being that broad, over many statistics,
#'     it is the least powerful of the families, and the targeted tests of
#'     \code{\link{ctScoreTest}} (stability, state dependence) find what they
#'     cover far more often. It costs a dual pass of the filter per eight
#'     parameters, and on the Laplace route one per truncated record, which is
#'     hours on large data.}
#'   \item{subjects, nonlinear, covariates, lagged}{Statistics of the data
#'     against datasets generated from the fit: how much subjects differ in
#'     mean, standard deviation and lag-one correlations; whether change
#'     depends on level nonlinearly; with \code{covariates}, whether subjects'
#'     statistics depend on them; and the data's own lagged correlations
#'     (\code{lags}). These use the data twice and are conservative; the
#'     lagged correlations especially, as the fit reproduces them, so they are
#'     not run unless asked for.}
#' }
#'
#' The error-based statistics are moment tests: each unit's contribution is
#' projected on its scores, which removes what estimating the parameters did
#' to the statistic (\code{projection}), and the variance is the errors' own:
#' they are independent over rows under the model, so a sum of products over
#' pairs of rows has the sum of the products' squares as its variance, which
#' holds with heavy tails or a changing variance; estimating the parameters
#' takes a known part of it away. That is the variance of the residual,
#' volatility and within-row families; the distribution and fit families,
#' whose statistics are sums of high powers, use the spread of the units'
#' contributions and a studentised bootstrap over units. Within a family the
#' largest standardised statistic and the sum of their squares are calibrated
#' with the statistics' correlation over units
#' and combined, with any exact p-values the family has, by the Cauchy
#' combination; families are combined by Holm. \code{adequate} is \code{TRUE}
#' when no family's adjusted p-value is below \code{alpha}.
#'
#' With many observations every misfit is significant, so the sizes of the
#' residual, latent and within-row correlations the model leaves are reported
#' too, as
#' root mean squares (the unbiased form of the SRMR, Maydeu-Olivares 2017),
#' with intervals and a test that the size is at most \code{close}.
#' \code{close_fit} is \code{TRUE} when every family that fails is one of those
#' three and its size is not significantly above \code{close}.
#'
#' Detection is not localization: a failure in one variable's equation reaches
#' the others' errors through the fitted covariance. \code{\link{ctScoreTest}}
#' is the more sensitive check for the departures it covers, and
#' \code{\link{ctModelImprove}} uses both.
#'
#' @param fit A fit from \code{\link{ctFit}}, either backend. The projection
#'   and the scores it needs require \code{backend = 'julia'}.
#' @param covariates Optional subject-level covariates, as for
#'   \code{\link{ctScoreTest}}, for the \code{covariates} family.
#' @param families Which families to run; see Details.
#' @param gaps Break points, in the data's time units, for the residual
#'   correlations' gap bins. By default the eighth, quarter and half quantiles
#'   of the time between every pair of a unit's rows.
#' @param close The root mean square residual correlation regarded as a close
#'   fit.
#' @param alpha Level for \code{adequate} and \code{close_fit}.
#' @param projection How the error-based statistics account for estimation:
#'   \code{'derivative'} differences the statistics along every parameter (two
#'   passes of the filter per parameter) and projects with the fit's Hessian;
#'   \code{'plugin'}, on a Laplace fit, takes those differences from plug-in
#'   passes, which filter each subject at effects solved once from the whole
#'   record rather than row by row -- several times cheaper, and approximate
#'   (on a pilot-shaped model it recovered the distribution family's misfit,
#'   z 5.3 against 5.6 exact and 1.8 uncorrected, and half the latent
#'   family's); \code{'opg'} uses the outer product of the scores, valid when
#'   units outnumber parameters many times; \code{'none'} ignores estimation,
#'   which is conservative, and for the distribution family can hide a strong
#'   misfit. \code{'auto'} takes the first that is affordable:
#'   \code{'derivative'} when the fit's Hessian is stored and the passes take
#'   under \code{getOption('ctsem.adequacy.budget', 120)} seconds, then
#'   \code{'plugin'} on the same terms, \code{'opg'} with at least ten units
#'   a parameter, else \code{'none'}.
#' @param parameters For the \code{scores} family, the raw parameters whose
#'   row scores are used: by default up to 16 of the dynamic and measurement
#'   parameters, those in DRIFT, DIFFUSION and CINT first.
#' @param nsamples Number of generated datasets, for the families that need
#'   them. Reused from \code{fit} when it already carries at least as many.
#' @param lags Observation lags for the \code{lagged} family.
#' @param nboot Bootstrap draws for each family's calibration.
#' @param cores Cores for the engine and for generating.
#' @return An object of class \code{ctFitAdequacy}: \code{families} (one row
#'   per family: the number of statistics, the worst, its z, the family's
#'   p-value and the Holm-adjusted one), \code{statistics} (every statistic:
#'   observed value, z, two-sided p), \code{size} (the root mean square
#'   residual correlations, their intervals and close-fit tests),
#'   \code{adequate}, \code{close_fit}, \code{projection} (the method used),
#'   and \code{fit} (with any generated data, for reuse).
#' @references Newey, W. K. (1985). Maximum likelihood specification testing
#'   and conditional moment tests. \emph{Econometrica, 53}, 1047-1070.
#'
#'   White, H. (1987). Specification testing in dynamic models. In T. Bewley
#'   (Ed.), \emph{Advances in Econometrics: Fifth World Congress}. Cambridge
#'   University Press.
#'
#'   Verbeke, G., & Molenberghs, G. (2013). The gradient function as an
#'   exploratory goodness-of-fit assessment of the random-effects distribution
#'   in mixed models. \emph{Biostatistics, 14}, 477-490.
#'
#'   Drikvandi, R., Verbeke, G., & Molenberghs, G. (2017). Diagnosing
#'   misspecification of the random-effects distribution in mixed models.
#'   \emph{Biometrics, 73}, 63-71.
#'
#'   Chang, J., Yao, Q., & Zhou, W. (2017). Testing for high-dimensional white
#'   noise using maximum cross-correlations. \emph{Biometrika, 104}, 111-127.
#'
#'   Liu, Y., & Xie, J. (2020). Cauchy combination test. \emph{Journal of the
#'   American Statistical Association, 115}, 393-402.
#'
#'   Maydeu-Olivares, A. (2017). Assessing the size of model misfit in
#'   structural equation models. \emph{Psychometrika, 82}, 533-558.
#' @seealso \code{\link{ctPostPredPlots}} and \code{\link{ctFitCheckCov}} to
#'   see the data against the fit.
#' @export
ctFitAdequacy <- function(fit, covariates = NULL,
  families = c('residuals', 'latent', 'volatility', 'pairs', 'distribution', 'fit',
    'dependence', 'effects', 'subjects', 'nonlinear', 'covariates'),
  gaps = NULL, close = 0.05, alpha = 0.05,
  projection = c('auto', 'derivative', 'plugin', 'opg', 'none'), parameters = NULL,
  nsamples = 100, lags = 0:3, nboot = 1999L, cores = 2) {
  if (!.ctFitIsJulia(fit) && !inherits(fit, 'ctStanFit')) stop('Not a ctsem fit.',
    call. = FALSE)
  asked <- !missing(families)
  families <- match.arg(families, c('residuals', 'latent', 'volatility', 'pairs',
    'distribution', 'fit', 'dependence', 'effects', 'scores', 'subjects', 'nonlinear',
    'covariates', 'lagged'), several.ok = TRUE)
  if (is.null(covariates)) families <- setdiff(families, 'covariates')
  # The effect distribution needs effects held as Laplace coordinates at one
  # level; by default it is run where that is so, and said where it is not.
  if ('effects' %in% families) {
    why <- .ctAdequacyEffectsUnavailable(fit)
    if (!is.null(why)) {
      if (asked) message("effects not tested: ", why)
      families <- setdiff(families, 'effects')
    }
  }
  if ('scores' %in% families && !.ctFitIsJulia(fit)) {
    message("scores not tested: row scores need backend = 'julia'.")
    families <- setdiff(families, 'scores')
  }
  projection <- match.arg(projection)
  if (.ctFitIsJulia(fit)) .ctBackendResolveThreads(cores)
  if (exists('.Random.seed', envir = globalenv(), inherits = FALSE)) {
    oldseed <- get('.Random.seed', envir = globalenv(), inherits = FALSE)
    on.exit(assign('.Random.seed', oldseed, envir = globalenv()), add = TRUE)
  } else on.exit(rm('.Random.seed', envir = globalenv()), add = TRUE)
  set.seed(20261010L)

  statistics <- list(); fams <- list(); size <- NULL; used <- NA_character_
  moments <- intersect(families, c('residuals', 'latent', 'volatility', 'pairs',
    'distribution', 'fit', 'effects', 'scores'))
  if (length(moments) || 'dependence' %in% families) {
    # The errors (on the Laplace route a sequential pass, the expensive part)
    # only for the families that read them: the effect distribution does not.
    neederr <- length(setdiff(c(moments, intersect(families, 'dependence')), 'effects')) > 0
    t0 <- proc.time()[['elapsed']]
    err <- if (neederr) .ctAdequacyErrors(fit) else NULL
    setup <- if (neederr) .ctAdequacySetup(err, gaps) else NULL
    if ('distribution' %in% moments && isTRUE(err$approximate)) message("The row ",
      "distances are not tested: on this fit each row's predictive covariance is an ",
      "approximation (random effects in the dynamics or variances, on the Laplace route).")
    if (length(moments)) {
      est <- if (.ctFitIsJulia(fit)) .ctFitRawEstimate(fit) else NULL
      extra <- list(
        directions = if ('scores' %in% moments) .ctAdequacyScoreDirections(fit, parameters) else NULL)
      # Every family's statistics at any parameter vector, so that the
      # projection can difference them all; at the estimate the errors are
      # already in hand.
      momentsAt <- function(x, errx = NULL, sequential = TRUE) {
        if (is.null(errx) && neederr) errx <- .ctAdequacyErrors(fit, x, sequential = sequential)
        .ctAdequacyMomentsAll(fit, x, errx, setup, moments, extra, cores)
      }
      mom <- momentsAt(est, err)
      tpass <- proc.time()[['elapsed']] - t0
      # Plug-in passes differ from the sequential ones only where the effects
      # are re-solved row by row: the Laplace route.
      laplace <- .ctFitIsJulia(fit) && identical(.ctBackendIntOverPop(fit$model_spec), 'laplace')
      proj <- .ctAdequacyProject(fit, mom, function(x) colSums(momentsAt(x)$M),
        projection, tpass, cores, plugin = if (laplace)
          function(x) colSums(momentsAt(x, sequential = FALSE)$M))
      used <- proj$method
      tests <- .ctAdequacyMomentTests(proj$Mstar, mom$V0, proj$reduction, mom$meta,
        nboot, raw = mom$M)
      extras <- .ctAdequacyExtras(err, moments)
      statistics[[length(statistics) + 1L]] <- data.frame(family = mom$meta$family,
        statistic = mom$meta$name, observed = tests$estimate, z = tests$z, p = tests$p,
        stringsAsFactors = FALSE)
      if (nrow(extras$statistics)) statistics[[length(statistics) + 1L]] <- extras$statistics
      for (f in moments) {
        r <- tests$families[[f]]
        if (is.null(r)) next
        ex <- extras$p[[f]]
        fams[[length(fams) + 1L]] <- data.frame(family = f,
          statistics = r$k + length(ex), worst = r$worst, worst_z = r$worst_z,
          p = .ctCauchyCombine(c(r$pmax, r$psum, ex)), stringsAsFactors = FALSE)
      }
      size <- .ctAdequacySize(tests, mom$meta, close)
    }
    if ('dependence' %in% families) {
      dep <- .ctAdequacyDependence(err)
      if (!is.null(dep)) {
        statistics[[length(statistics) + 1L]] <- dep$statistics
        fams[[length(fams) + 1L]] <- dep$family
      }
    }
  }

  generated <- intersect(families, c('subjects', 'nonlinear', 'covariates', 'lagged'))
  if (length(generated)) {
    gen <- .ctAdequacyGenerated(fit, covariates, generated, nsamples, lags, cores)
    fit <- gen$fit
    statistics[[length(statistics) + 1L]] <- gen$statistics
    fams <- c(fams, list(gen$families))
  }
  statistics <- do.call(rbind, statistics)
  families_out <- do.call(rbind, fams)
  rownames(statistics) <- rownames(families_out) <- NULL
  families_out$p_adj <- stats::p.adjust(families_out$p, 'holm')
  failing <- families_out$family[families_out$p_adj < alpha]
  sized <- if (is.null(size)) character() else
    size$family[!is.na(size$p_close) & size$p_close >= alpha & size$overall]
  structure(list(families = families_out, statistics = statistics, size = size,
    adequate = !length(failing), close_fit = all(failing %in% sized),
    alpha = alpha, close = close, projection = used,
    nsamples = if (length(generated)) nsamples else NA_integer_, fit = fit),
    class = 'ctFitAdequacy')
}

# The one-step-ahead errors at a raw parameter vector (the estimate unless
# given): each variable's error given the row's other observed errors,
# standardised -- (P v)_j / sqrt(P_jj), P the inverse of the observed block of
# the predicted covariance -- the row's squared whitened distance v' P v, its
# dimension, and the correlation the model implies between two variables'
# conditional errors in the row: P v has covariance P, so it is
# P_jk / sqrt(P_jj P_kk) -- the partial correlation of the observations with
# its sign reversed, which a first version used and which rejected every
# correctly specified model. On the Laplace route
# the predictions are each unit's from its earlier rows alone (sequential), as
# the errors must be for their independence over rows to hold.
.ctAdequacyErrors <- function(fit, raw = NULL, sequential = TRUE) {
  model <- .ctFitModelObject(fit)
  vars <- model$manifestNames
  m <- length(vars)
  if (.ctFitIsJulia(fit)) {
    if (is.null(raw)) raw <- .ctFitRawEstimate(fit)
    spec <- .ctBackendSpec(fit)
    k <- suppressMessages(.ctBackendKalmanRaw(fit, raw, subjectmatrices = FALSE,
      fields = c('eta', 'etacov', 'y', 'ycov', 'subject'), sequential = sequential))
    nrows <- length(spec$times)
    Y <- matrix(t(spec$manifest_data), nrows, m)
    yp <- matrix(k$y[1, , , drop = FALSE], nrows, m)
    yc <- array(k$ycov[1, , , , drop = FALSE], c(nrows, m, m))
    naug <- dim(k$eta)[3]
    nl <- spec$nlatent
    eta <- matrix(k$eta[1, , , drop = FALSE], nrows, naug)[, seq_len(nl), drop = FALSE]
    etaupd <- matrix(k$eta[2, , , drop = FALSE], nrows, naug)[, seq_len(nl), drop = FALSE]
    covdiff <- vapply(seq_len(nl), function(j) k$etacov[1, , j, j] - k$etacov[2, , j, j],
      numeric(nrows))
    unit <- as.integer(k$subject)
    time <- as.numeric(spec$times)
  } else {
    if (!is.null(raw)) stop('Stan fits are tested at their estimate only.', call. = FALSE)
    k <- suppressMessages(ctKalmanArray(fit))
    nrows <- nrow(k$y)
    Y <- matrix(k$y, nrows, m)
    yp <- matrix(k$yprior[1, , , drop = FALSE], nrows, m)
    yc <- array(k$ypriorcov[1, , , , drop = FALSE], c(nrows, m, m))
    eta <- matrix(k$etaprior[1, , , drop = FALSE], nrows)
    nl <- ncol(eta)
    etaupd <- matrix(k$etaupd[1, , , drop = FALSE], nrows)[, seq_len(nl), drop = FALSE]
    covdiff <- vapply(seq_len(nl), function(j) k$etapriorcov[1, , j, j] -
      k$etaupdcov[1, , j, j], numeric(nrows))
    unit <- as.integer(factor(k$id, levels = unique(k$id)))
    time <- as.numeric(k$time)
  }
  Y[!is.finite(Y) | Y == 99999] <- NA
  E <- matrix(NA_real_, nrows, m, dimnames = list(NULL, vars))
  R <- array(NA_real_, c(nrows, m, m))
  Q <- rep(NA_real_, nrows); d <- integer(nrows)
  for (r in seq_len(nrows)) {
    obs <- which(!is.na(Y[r, ]))
    if (!length(obs)) next
    v <- Y[r, obs] - yp[r, obs]
    P <- try(solve(matrix(yc[r, obs, obs], length(obs))), silent = TRUE)
    if (inherits(P, 'try-error') || any(!is.finite(P))) next
    Pv <- as.numeric(P %*% v)
    s <- sqrt(diag(P))
    E[r, obs] <- Pv / s
    R[r, obs, obs] <- P / outer(s, s)
    Q[r] <- sum(v * Pv); d[r] <- length(obs)
  }
  # Each latent's update by the row, standardised: the state moved by K v,
  # whose covariance under the model is the prior less the updated -- a
  # function of the row's innovation alone given the past, so uncorrelated
  # across rows like the errors, but one per latent process rather than one
  # per item, which is where a structural misfit shows.
  covdiff <- matrix(covdiff, nrows)
  Z <- (etaupd - eta) / sqrt(pmax(covdiff, 0))
  Z[!is.finite(Z) | covdiff <= 1e-12 * pmax(1, abs(covdiff))] <- NA
  colnames(Z) <- model$latentNames[seq_len(ncol(Z))]
  labels <- tryCatch(as.character(.ctFitIdMap(fit)$original), error = function(e) NULL)
  if (length(labels) < max(unit)) labels <- as.character(seq_len(max(unit)))
  list(E = E, R = R, Q = Q, d = d, Z = Z, unit = unit, nunits = max(unit), time = time,
    eta = eta, vars = vars, labels = labels, approximate = .ctAdequacyApproximate(fit))
}

# Whether each row's predictive covariance is a Gaussian approximation: on the
# Laplace route, where the subject's random effects are integrated by their
# posterior given the earlier rows, the sequential predictions add that
# posterior's uncertainty by the delta method, which is exact only where the
# effects enter the predictions linearly -- intercepts, initial means,
# predictor effects. Where they reach the dynamics or the variances, the true
# predictive is a mixture, a little wider. Generated from a pilot-shaped model
# (random drift and diffusion, 4500 rows), the row distances ran 1.5 per cent
# high, and with the estimation removed their mean rejected the correct model
# at z = 6.0 (z = 1.7 without removing it).
.ctAdequacyApproximate <- function(fit) {
  if (!.ctFitIsJulia(fit) || !identical(.ctBackendIntOverPop(fit$model_spec), 'laplace'))
    return(FALSE)
  pars <- .ctFitModelObject(fit)$pars
  varying <- pars$matrix[as.logical(pars$indvarying) %in% TRUE & !is.na(pars$param)]
  any(!varying %in% c('MANIFESTMEANS', 'CINT', 'T0MEANS', 'TDPREDEFFECT'))
}

# What stays fixed when the statistics are recomputed at other parameter
# values: the gap bins, and the mean the squared errors are centred on.
# Centring each variable's squares on their overall mean leaves the volatility
# statistic about clustering in time -- a variance the model gets wrong
# everywhere is the distribution family's -- and costs no estimation effect,
# as the statistic's derivative in that mean is zero under the model.
.ctAdequacySetup <- function(err, gaps = NULL) {
  if (is.null(gaps)) {
    rows <- split(seq_along(err$unit), err$unit)
    g <- unlist(lapply(rows, function(r) {
      if (length(r) < 2L) return(numeric())
      t <- err$time[r]
      o <- outer(t, t, function(a, b) b - a)
      o[upper.tri(o)]
    }), use.names = FALSE)
    gaps <- if (length(g) > 10L) unique(stats::quantile(g, c(1 / 8, 1 / 4, 1 / 2),
      names = FALSE)) else numeric()
    gaps <- gaps[gaps > 0]
  }
  gaps <- sort(unique(as.numeric(gaps)))
  edges <- c(-Inf, gaps, Inf)
  fmt <- function(x) format(signif(x, 3))
  labels <- if (!length(gaps)) 'all' else c(paste0('<=', fmt(gaps[1])),
    if (length(gaps) > 1L) paste0(fmt(gaps[-length(gaps)]), '-', fmt(gaps[-1])),
    paste0('>', fmt(gaps[length(gaps)])))
  list(edges = edges, labels = labels,
    centre = colMeans(err$E^2, na.rm = TRUE))
}

# Each unit's contribution to every moment statistic, a units x statistics
# matrix, and what each column is. Missing errors contribute nothing, which
# leaves every contribution's mean at zero under the model.
#
# `V0` holds each unit's contribution to the statistic's variance under the
# model, from the rows and row pairs themselves: the errors are independent
# over rows, so a sum of products over pairs of rows has variance the sum of
# the products' squares (the form Escanciano & Lobato 2009 use, robust to
# heavy tails and to variance that changes over time), and a sum over rows the
# sum of the terms' squares. The spread of the units' totals estimates the
# same thing under the model, but it also allows the dependence within a unit
# that the null rules out, so under the alternative it grows with exactly the
# clustering being tested: a volatility family on it found state-dependent
# diffusion 16 times in 100 where the pairwise variance finds it 88. The
# unit spread is kept where the unit is the statistic's grain (`fit`, NA here)
# and for the statistics' correlation.
.ctAdequacyMoments <- function(err, setup, which) {
  E <- err$E; p <- ncol(E); n <- err$nunits; vars <- err$vars
  unitsum <- function(v) {
    out <- matrix(0, n, NCOL(v))
    s <- rowsum(as.matrix(v), err$unit, reorder = TRUE)
    out[as.integer(rownames(s)), ] <- s
    out
  }
  blocks <- list()
  addblock <- function(M, V0, family, name, count, kind, bin = NA_integer_)
    blocks[[length(blocks) + 1L]] <<- list(M = M, V0 = V0, meta = data.frame(
      family = family, name = name, count = count, kind = kind, bin = bin,
      stringsAsFactors = FALSE))
  if (any(c('residuals', 'volatility') %in% which)) {
    B <- length(setup$labels)
    if ('residuals' %in% which) {
      blk <- .ctAdequacyPairBlock(E, err$unit, err$time, n, setup)
      grid <- expand.grid(j = seq_len(p), k = seq_len(p), b = seq_len(B))
      addblock(blk$M, blk$V0, 'residuals', paste0(vars[grid$j], '~', vars[grid$k], ' gap ',
        setup$labels[grid$b]), blk$N, 'corr', grid$b)
    }
    if ('volatility' %in% which) {
      E2 <- sweep(E^2, 2, setup$centre)
      blk <- .ctAdequacyPairBlock(E2, err$unit, err$time, n, setup)
      grid <- expand.grid(j = seq_len(p), k = seq_len(p), b = seq_len(B))
      own <- grid$j == grid$k
      addblock(blk$M[, own, drop = FALSE], blk$V0[, own, drop = FALSE], 'volatility',
        paste0(vars[grid$j[own]], '^2 gap ', setup$labels[grid$b[own]]), blk$N[own],
        'square', grid$b[own])
    }
  }
  if ('latent' %in% which) {
    q <- ncol(err$Z)
    blk <- .ctAdequacyPairBlock(err$Z, err$unit, err$time, n, setup)
    grid <- expand.grid(j = seq_len(q), k = seq_len(q), b = seq_len(length(setup$labels)))
    lat <- colnames(err$Z)
    addblock(blk$M, blk$V0, 'latent', paste0(lat[grid$j], '~', lat[grid$k], ' gap ',
      setup$labels[grid$b]), blk$N, 'corr', grid$b)
  }
  if ('pairs' %in% which && p > 1L) {
    pr <- which(upper.tri(diag(p)), arr.ind = TRUE)
    pr <- pr[order(pr[, 1], pr[, 2]), , drop = FALSE]
    V <- sapply(seq_len(nrow(pr)), function(a) {
      v <- E[, pr[a, 1]] * E[, pr[a, 2]] - err$R[, pr[a, 1], pr[a, 2]]
      v[!is.finite(v)] <- 0
      v
    })
    V <- matrix(V, nrow(E))
    N <- colSums(!is.na(E[, pr[, 1], drop = FALSE]) & !is.na(E[, pr[, 2], drop = FALSE]))
    addblock(unitsum(V), unitsum(V^2), 'pairs', paste0(vars[pr[, 1]], ':', vars[pr[, 2]]),
      N, 'corr')
  }
  if ('distribution' %in% which) {
    E3 <- E^3; E3[is.na(E3)] <- 0
    E4 <- E^4 - 3; E4[is.na(E4)] <- 0
    q1 <- err$Q - err$d; q1[!is.finite(q1)] <- 0
    q2 <- (err$Q - err$d)^2 - 2 * err$d; q2[!is.finite(q2)] <- 0
    nobs <- colSums(!is.na(E))
    # The row distances test the predictive covariance's level and spread,
    # which on an approximate route (.ctAdequacyApproximate) they find wrong
    # in the approximation rather than the model.
    distances <- !isTRUE(err$approximate)
    terms <- if (distances) cbind(E3, E4, q1, q2) else cbind(E3, E4)
    # Unit-level variance and the studentised bootstrap: a sum of fourth powers
    # has a pairwise-style variance in eighth powers, too noisy for a normal
    # reference (4.7 per cent of correct models rejected at .01), and rows stay
    # independent under the departures this family is for, so nothing is lost.
    addblock(unitsum(terms), matrix(NA_real_, n, ncol(terms)), 'distribution',
      c(paste0('skew(', vars, ')'), paste0('kurtosis(', vars, ')'),
        if (distances) c('distance', 'distance spread')),
      c(nobs, nobs, if (distances) rep(sum(err$d > 0), 2)), 'moment')
  }
  if ('fit' %in% which) {
    Q <- err$Q; Q[!is.finite(Q)] <- 0
    qi <- unitsum(Q); di <- unitsum(err$d)
    addblock((qi - di)^2 - 2 * di, matrix(NA_real_, n, 1), 'fit', 'spread of unit fit',
      n, 'moment')
  }
  list(M = do.call(cbind, lapply(blocks, `[[`, 'M')),
    V0 = do.call(cbind, lapply(blocks, `[[`, 'V0')),
    meta = do.call(rbind, lapply(blocks, `[[`, 'meta')))
}

# Sums over each unit's pairs of rows s < t, binned by the time between them,
# of products of the columns of X (rows by K, in the units' row order): the
# statistic for every column pair and bin, per unit (`M`, units by K*K*B, the
# first index fastest), its variance under the model from the products'
# squares (`V0`), and the number of pairs that contribute (`N`). Missing
# values contribute nothing.
.ctAdequacyPairBlock <- function(X, unit, time, n, setup) {
  X <- as.matrix(X)
  K <- ncol(X)
  B <- length(setup$labels)
  Xz <- X; Xz[!is.finite(Xz)] <- 0
  O <- is.finite(X) * 1
  X2 <- Xz^2
  M <- V0 <- matrix(0, n, K * K * B); N <- numeric(K * K * B)
  rows <- split(seq_along(unit), factor(unit, levels = seq_len(n)))
  for (i in seq_len(n)) {
    r <- rows[[i]]
    if (length(r) < 2L) next
    t <- time[r]
    bin <- matrix(.bincode(outer(t, t, function(a, b) b - a), setup$edges, right = TRUE,
      include.lowest = TRUE), length(r))
    bin[lower.tri(bin, diag = TRUE)] <- 0L
    x <- Xz[r, , drop = FALSE]; x2 <- X2[r, , drop = FALSE]; o <- O[r, , drop = FALSE]
    for (b in seq_len(B)) {
      W <- (bin == b) * 1
      if (!any(W > 0)) next
      cols <- (b - 1L) * K * K + seq_len(K * K)
      M[i, cols] <- as.numeric(crossprod(x, W %*% x))
      V0[i, cols] <- as.numeric(crossprod(x2, W %*% x2))
      N[cols] <- N[cols] + as.numeric(crossprod(o, W %*% o))
    }
  }
  list(M = M, V0 = V0, N = N)
}

# Every requested family's per-unit statistics at the parameters `x` (NULL:
# the estimate of a fit that has no other), from the errors `err` there.
.ctAdequacyMomentsAll <- function(fit, x, err, setup, which, extra, cores) {
  parts <- list()
  errfams <- intersect(which, c('residuals', 'latent', 'volatility', 'pairs', 'distribution',
    'fit'))
  if (length(errfams)) parts[[length(parts) + 1L]] <- .ctAdequacyMoments(err, setup, errfams)
  if ('scores' %in% which) parts[[length(parts) + 1L]] <- .ctAdequacyScoreMoments(fit, x,
    err, setup, extra$directions, cores)
  if ('effects' %in% which) parts[[length(parts) + 1L]] <- .ctAdequacyEffectMoments(fit,
    x, cores)
  list(M = do.call(cbind, lapply(parts, `[[`, 'M')),
    V0 = do.call(cbind, lapply(parts, `[[`, 'V0')),
    meta = do.call(rbind, lapply(parts, `[[`, 'meta')))
}

# Why the effect distribution cannot be tested on this fit, or NULL.
.ctAdequacyEffectsUnavailable <- function(fit) {
  if (!.ctFitIsJulia(fit)) return("it needs backend = 'julia'.")
  if (!identical(.ctBackendIntOverPop(fit$model_spec), 'laplace'))
    return("it needs random effects integrated by intoverpop = 'laplace'.")
  dims <- try(.ctBackendUnitDims(fit), silent = TRUE)
  if (inherits(dims, 'try-error')) return("the effects are not at a single level.")
  if (!length(dims) || all(dims == 0)) return("the model has no random effects.")
  if (length(unique(dims)) > 1L) return("the units differ in their number of effects.")
  if (dims[1] > 4L) return(paste0(dims[1], " effects a unit; the posterior by ",
    "quadrature would take ", 4^dims[1], " evaluations a unit."))
  NULL
}

# The effects' distribution, unit by unit: the posterior expectation of
# Phi((c - u_k) / h) along each axis k of the unit's standardised effects, less
# the prior's, Phi(c / sqrt(1 + h^2)) -- a smoothed posterior probability below
# each cut. Averaged over units each has mean zero under the model, since the
# posterior averages to the prior, and each is bounded. The gradient function
# these integrate, f(y | u) / m(y) at a point u, is not: in the tail it is near
# zero for almost every unit and huge for the rare one near the point, and a
# studentised mean of it rejected a correct model a third of the time at .05.
# The posterior is the quadrature's own rule, whose error for these at h = 0.75
# is at most 1e-4 a unit for one effect (nine nodes) and 3e-3 for three (five),
# even with a posterior as wide as the prior.
.ctAdequacyEffectMoments <- function(fit, x, cores) {
  if (is.null(x)) x <- .ctFitRawEstimate(fit)
  d <- .ctBackendUnitDims(fit)[1]
  post <- .ctBackendPosteriorNodes(fit, x, nodes = c(9L, 7L, 5L, 4L)[min(d, 4L)], cores)
  W <- post$weights
  h <- 0.75
  cuts <- c(-1.5, -0.75, 0, 0.75, 1.5)
  C <- do.call(cbind, lapply(seq_len(d), function(k) {
    Pk <- matrix(post$points[, , k], nrow(W))
    vapply(cuts, function(cc) rowSums(W * stats::pnorm((cc - Pk) / h)) -
      stats::pnorm(cc / sqrt(1 + h^2)), numeric(nrow(W)))
  }))
  C[!is.finite(C)] <- 0
  n <- nrow(C)
  list(M = C, V0 = matrix(NA_real_, n, ncol(C)), meta = data.frame(family = 'effects',
    name = paste0('P(u', rep(seq_len(d), each = length(cuts)), ' < ', cuts, ')'),
    count = n, kind = 'moment', bin = NA_integer_, stringsAsFactors = FALSE))
}

# The raw parameters whose row scores the scores family uses: named, or up to
# 16 of the fit's dynamic and measurement parameters, DRIFT, DIFFUSION and
# CINT first. As unit directions, columns named.
.ctAdequacyScoreDirections <- function(fit, parameters = NULL) {
  ptab <- .ctScoreParameterTable(fit)
  keep <- ptab$kind == 'model' & !ptab$matrix %in% c('T0MEANS', 'T0VAR')
  if (!is.null(parameters)) {
    unknown <- setdiff(parameters, ptab$name)
    if (length(unknown)) stop("Unknown parameters: ", paste(unknown, collapse = ', '),
      call. = FALSE)
    keep <- ptab$name %in% parameters
  }
  idx <- which(keep)
  first <- ptab$matrix[idx] %in% c('DRIFT', 'DIFFUSION', 'CINT')
  idx <- utils::head(c(idx[first], idx[!first]), if (is.null(parameters)) 16L else length(idx))
  D <- diag(nrow(ptab))[, idx, drop = FALSE]
  colnames(D) <- ptab$name[idx]
  D
}

# White's dynamic information matrix test: each row's score for the chosen
# parameters, multiplied by the scores of the same unit's later rows, by
# time gap. Under the model a unit's row scores are a martingale difference
# sequence, so every such product has mean zero.
#
# Unit-level, like the distribution family, rather than with the pairwise
# variance the errors use. A row score's variance given the past is not
# constant -- a drift score scales with the predicted state -- and the per-pair
# variance, which assumes it is, rejected a correct model a quarter of the time
# at .05. The form that needs only the martingale property, each later row's
# squared score times its earlier partners' summed score squared, got the
# median variance right (ratio to the projected spread 0.96-1.05 over six null
# fits) but not the tails: products of chi-square-like variance scores over
# long gaps are far from normal, and where estimation takes half a statistic's
# variance the difference of two noisy terms is fragile, so the largest of 64
# |z| reached 3.5-4.7. The units' spread and the studentised bootstrap over
# units see both.
.ctAdequacyScoreMoments <- function(fit, x, err, setup, directions, cores) {
  if (is.null(x)) x <- .ctFitRawEstimate(fit)
  R <- .ctBackendRowScores(fit, x, directions, cores)
  blk <- .ctAdequacyPairBlock(R, err$unit, err$time, err$nunits, setup)
  K <- ncol(R)
  grid <- expand.grid(j = seq_len(K), k = seq_len(K), b = seq_len(length(setup$labels)))
  par <- colnames(directions)
  list(M = blk$M, V0 = matrix(NA_real_, nrow(blk$M), ncol(blk$M)),
    meta = data.frame(family = 'scores',
    name = paste0('s(', par[grid$j], ')~s(', par[grid$k], ') gap ', setup$labels[grid$b]),
    count = blk$N, kind = 'moment', bin = grid$b, stringsAsFactors = FALSE))
}

# Each unit's contribution with the estimation effect removed:
# m_i + Gamma I^-1 s_i, Gamma the derivative of the summed statistics in the
# parameters, I the information, s_i the unit's scores (Newey 1985). Summed,
# the contributions also correct for an estimate that is not exactly a root of
# the scores, as a quadrature-continued Laplace estimate is not. Estimating
# the parameters takes Gamma I^-1 Gamma' out of each statistic's variance under
# the model (`reduction`; Box & Pierce's correction for residual
# autocorrelations is the special case).
#
# Gamma by central differences of the statistics along every parameter is
# exact, and costs two filter passes a parameter -- seconds on a small model,
# hours on the SNSF pilot, where one sequential Laplace pass took 30 minutes.
# The outer-product form, Gamma = -sum m_i s_i' with I = sum s_i s_i', needs
# nothing more than the scores, but regresses each statistic on the p scores
# over n units, which with p near n removes variance that is not the
# estimation's; it is used only at ten units a parameter.
.ctAdequacyProject <- function(fit, mom, sums, projection, tpass, cores, plugin = NULL) {
  n <- nrow(mom$M)
  none <- function(note = NULL) {
    if (!is.null(note)) message(note)
    list(Mstar = mom$M, reduction = numeric(ncol(mom$M)), method = 'none')
  }
  if (identical(projection, 'none')) return(none())
  if (!.ctFitIsJulia(fit)) return(none(if (projection != 'auto')
    "Stan fits are tested without projection; the references are conservative."))
  est <- .ctFitRawEstimate(fit)
  p <- length(est)
  S <- try(.ctBackendWithMaxChunks(cores, .ctBackendScoreMatrix(fit, est)), silent = TRUE)
  if (inherits(S, 'try-error') || nrow(S) != n) return(none(paste0("The units' ",
    "scores are not available for this model; the error-based tests ignore ",
    "estimation and are conservative.")))
  budget <- getOption('ctsem.adequacy.budget', 120)
  stored <- !is.null(.ctBackendStoredHessian(fit, est))
  method <- projection
  if (method == 'plugin' && is.null(plugin)) method <- 'derivative'
  # Gamma by central differences of `f` along every parameter, or NULL when
  # `auto` finds the first parameter's pair would put the whole over budget --
  # timed there rather than on the base pass, which pays the compile on a
  # session's first call.
  differences <- function(f) {
    G <- matrix(0, ncol(mom$M), p)
    for (k in seq_len(p)) {
      started <- proc.time()[['elapsed']]
      h <- 1e-3 * max(1, abs(est[k]))
      e <- est; e[k] <- est[k] + h
      plus <- f(e)
      e[k] <- est[k] - h
      minus <- f(e)
      G[, k] <- (plus - minus) / (2 * h)
      if (k == 1L && method == 'auto' &&
          (proc.time()[['elapsed']] - started) * p > budget) return(NULL)
    }
    G
  }
  Gamma <- NULL; used <- NA_character_
  if (method %in% c('auto', 'derivative') && (method == 'derivative' || stored)) {
    Gamma <- differences(sums); used <- 'derivative'
  }
  if (is.null(Gamma) && !is.null(plugin) && (method == 'plugin' ||
      (method == 'auto' && stored))) {
    Gamma <- differences(plugin); used <- 'plugin'
  }
  if (!is.null(Gamma)) {
    H <- .ctBackendWithMaxChunks(cores, .ctBackendHessian(fit, est))
    A <- Gamma %*% .ctScorePinv(-(H + t(H)) / 2)$inv
    return(list(Mstar = mom$M + S %*% t(A), reduction = rowSums(A * Gamma),
      method = used))
  }
  if (method %in% c('auto', 'opg') && (method == 'opg' || n >= 10L * p)) {
    coef <- .ctScorePinv(crossprod(S))$inv %*% crossprod(S, mom$M)
    fitted <- S %*% coef
    return(list(Mstar = mom$M - fitted, reduction = colSums(fitted * mom$M),
      method = 'opg'))
  }
  none(paste0("The error-based tests ignore estimation (", p, " parameters, ", n,
    " units; projection = 'derivative' would take about ", ceiling(2 * p * tpass / 60),
    " minutes", if (!is.null(plugin)) paste0(", 'plugin' several times less, and is ",
      "approximate") else "", ") and are conservative; the distribution family can miss ",
    "a strong misfit without it."))
}

# Every statistic's z, and each family's max and sum tests. A statistic's
# variance is its pairwise or row-wise variance under the model less what
# estimation took (`V0`, `reduction`), or for a unit-level statistic the
# spread of the units' projected contributions. A family's largest |z| and sum
# of squares are referred to draws of the standardised statistics with their
# correlation, taken from the units' projected contributions: Gaussian draws
# where the variances are pairwise or row-wise, sums over thousands of nearly
# independent terms; for a unit-level statistic, which is a skewed sum over a
# few hundred units, a studentised bootstrap over units, which corrects the
# skewness to second order (a Gaussian one rejected a correct model 14 per
# cent of the time at .05).
.ctAdequacyMomentTests <- function(Mstar, V0, reduction, meta, nboot = 1999L,
  raw = Mstar) {
  n <- nrow(Mstar)
  U <- colSums(Mstar)
  C <- sweep(Mstar, 2, colMeans(Mstar))
  spread <- colSums(C^2) * n / (n - 1)
  v0 <- colSums(V0)
  rowlevel <- is.finite(v0)
  v <- ifelse(rowlevel, v0 - reduction, spread)
  # A statistic the scores already determine -- the summed row distance is
  # nearly the score of the variances' common scale -- keeps almost none of
  # its variance once estimation is taken out, and a near-zero sum over a
  # near-zero variance gave z = 10 on correctly specified models. Under one
  # per cent left, it says nothing the fit did not, and is dropped.
  rawspread <- colSums(sweep(raw, 2, colMeans(raw))^2) * n / (n - 1)
  base <- ifelse(rowlevel, v0, rawspread)
  ok <- is.finite(v) & is.finite(base) & base > 0 & v > 0.01 * base
  z <- ifelse(ok, U / sqrt(pmax(v, 0)), NA_real_)
  p <- 2 * stats::pnorm(-abs(z))
  okcols <- which(ok)
  Ub <- matrix(NA_real_, nboot, length(okcols))
  families <- lapply(stats::setNames(nm = unique(meta$family)), function(f) {
    idx <- okcols[meta$family[okcols] == f]
    if (!length(idx)) return(NULL)
    at <- match(idx, okcols)
    zf <- z[idx]
    if (all(rowlevel[idx])) {
      R <- stats::cov2cor(crossprod(C[, idx, drop = FALSE]) + diag(1e-12, length(idx)))
      e <- eigen(R, symmetric = TRUE)
      keep <- e$values > 1e-8 * e$values[1]
      Zb <- matrix(stats::rnorm(nboot * sum(keep)), nboot) %*%
        (sqrt(e$values[keep]) * t(e$vectors[, keep, drop = FALSE]))
    } else {
      Zb <- matrix(NA_real_, nboot, length(idx))
      Cf <- C[, idx, drop = FALSE]
      for (b in seq_len(nboot)) {
        X <- Cf[sample.int(n, n, replace = TRUE), , drop = FALSE]
        u <- colSums(X)
        vb <- (colSums(X^2) - u^2 / n) * n / (n - 1)
        Zb[b, ] <- ifelse(vb > 0, u / sqrt(vb), NA_real_)
      }
    }
    Ub[, at] <<- sweep(Zb, 2, sqrt(v[idx]), '*')
    tmax <- max(abs(zf)); tsum <- sum(zf^2)
    bmax <- apply(abs(Zb), 1, max, na.rm = TRUE)
    bsum <- rowSums(Zb^2, na.rm = TRUE)
    w <- idx[which.max(abs(zf))]
    list(k = length(idx), worst = meta$name[w], worst_z = z[w],
      pmax = (1 + sum(bmax >= tmax)) / (nboot + 1),
      psum = (1 + sum(bsum >= tsum)) / (nboot + 1))
  })
  # The estimate each statistic is about: a correlation for the residual and
  # within-row statistics (the projected sum over the pairs that contribute),
  # the unit mean otherwise.
  denom <- ifelse(meta$kind %in% c('corr', 'square'), pmax(meta$count, 1), n)
  list(z = z, p = p, estimate = U / denom, U = U, v = v, ok = ok, denom = denom,
    families = families, draws = Ub)
}

# P-values the Cauchy combination takes beside the bootstrap's: the
# Anderson-Darling tests of the errors against the normal and of the row
# distances against the chi-square, and the worst unit's distance. Their
# references ignore estimation, which only makes them conservative.
.ctAdequacyExtras <- function(err, which) {
  stats <- list(); p <- list()
  zp <- function(x) stats::qnorm(pmax(x, 1e-300) / 2, lower.tail = FALSE)
  add <- function(family, name, observed, pv) {
    stats[[length(stats) + 1L]] <<- data.frame(family = family, statistic = name,
      observed = observed, z = zp(pv), p = pv, stringsAsFactors = FALSE)
    p[[family]] <<- c(p[[family]], pv)
  }
  if ('distribution' %in% which) {
    for (j in seq_along(err$vars)) {
      ad <- .ctAdequacyAD(stats::pnorm(err$E[, j]))
      if (!is.null(ad)) add('distribution', paste0('AD normal(', err$vars[j], ')'), ad[1], ad[2])
    }
    ok <- err$d > 0 & is.finite(err$Q)
    ad <- .ctAdequacyAD(stats::pchisq(err$Q[ok], err$d[ok]))
    if (!is.null(ad)) add('distribution', 'AD chisq(row distance)', ad[1], ad[2])
  }
  if ('fit' %in% which) {
    Q <- err$Q; Q[!is.finite(Q)] <- 0
    qi <- tapply(Q, err$unit, sum); di <- tapply(err$d, err$unit, sum)
    ok <- di > 0
    u <- stats::pchisq(qi[ok], di[ok])
    ext <- 2 * pmin(u, 1 - u)
    w <- which.min(ext)
    if (length(w)) add('fit', paste0('worst unit: ', err$labels[as.integer(names(u)[w])]),
      qi[ok][w] / di[ok][w],
      1 - (1 - ext[w])^length(u))
  }
  list(statistics = if (length(stats)) do.call(rbind, stats) else
    data.frame(family = character(), statistic = character(), observed = numeric(),
      z = numeric(), p = numeric()), p = p)
}

# The size of what the model leaves in the error correlations: the root mean
# square residual correlation, unbiased (each squared correlation less its
# sampling variance, Maydeu-Olivares 2017), with a delta-method interval and a
# one-sided test that it is at most `close`, the variance of the mean square
# from the bootstrap's draws of the summed contributions. For the residual correlations over all gaps, over
# the shortest gaps alone, and for the within-row correlations.
.ctAdequacySize <- function(tests, meta, close) {
  okcols <- which(tests$ok)
  measure <- function(label, family, cols) {
    cols <- intersect(cols, okcols)
    if (!length(cols)) return(NULL)
    r <- tests$U[cols] / tests$denom[cols]
    se <- sqrt(tests$v[cols]) / tests$denom[cols]
    k <- length(cols)
    ms <- mean(r^2 - se^2)
    # d(ms)/dU_k = 2 r_k / (k denom_k); the draws are U's null spread.
    at <- match(cols, okcols)
    lin <- tests$draws[, at, drop = FALSE] %*% (2 * r / (k * tests$denom[cols]))
    sems <- stats::sd(lin)
    est <- sqrt(max(ms, 0))
    data.frame(family = family, measure = label, estimate = est,
      lower = sqrt(max(ms - 1.96 * sems, 0)), upper = sqrt(max(ms + 1.96 * sems, 0)),
      close = close, p_close = stats::pnorm((ms - close^2) / sems, lower.tail = FALSE),
      overall = !grepl('shortest', label), stringsAsFactors = FALSE)
  }
  res <- which(meta$family == 'residuals' & meta$kind == 'corr')
  lat <- which(meta$family == 'latent')
  out <- list(
    measure('residual correlation, all gaps', 'residuals', res),
    measure('residual correlation, shortest gaps', 'residuals',
      res[meta$bin[res] %in% 1L]),
    measure('latent residual correlation, all gaps', 'latent', lat),
    measure('latent residual correlation, shortest gaps', 'latent',
      lat[meta$bin[lat] %in% 1L]),
    measure('within-row correlation', 'pairs', which(meta$family == 'pairs')))
  out <- do.call(rbind, out)
  if (!is.null(out)) rownames(out) <- NULL
  out
}

# p-values combined by the Cauchy combination (Liu & Xie 2020): the mean of
# tan((0.5 - p) pi) is standard Cauchy in the tail under any dependence.
.ctCauchyCombine <- function(p) {
  p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  t <- mean(tan((0.5 - p) * pi))
  0.5 - atan(t) / pi
}

# Distance correlation of each variable's error, and of its square, with what
# came before -- the predicted states, the previous row's errors, the time
# since it -- by permutation over rows. Under the model each error is
# independent of all of that, so the permutation reference is exact for known
# parameters (checked on exact Kalman innovations: 3.4 per cent at .05).
.ctAdequacyDependence <- function(err, nperm = 299L, nmax = 800L) {
  E <- err$E; vars <- err$vars; p <- length(vars); n <- nrow(E)
  unit <- err$unit
  prev <- c(NA_integer_, seq_len(n - 1L))
  prev[c(TRUE, unit[-1] != unit[-n])] <- NA_integer_
  X <- cbind(err$eta, E[prev, , drop = FALSE], dt = err$time - err$time[pmax(prev, 1L)])
  use <- which(!is.na(prev) & stats::complete.cases(X) & stats::complete.cases(E))
  if (length(use) < 50L) return(NULL)
  if (length(use) > nmax) use <- sort(sample(use, nmax))
  Xs <- scale(X[use, , drop = FALSE])
  Xs <- Xs[, apply(Xs, 2, function(x) all(is.finite(x))), drop = FALSE]
  A <- .ctAdequacyCentredDistance(Xs)
  targets <- list(); tn <- character()
  for (j in seq_len(p)) {
    targets[[length(targets) + 1L]] <- .ctAdequacyCentredDistance(E[use, j, drop = FALSE])
    tn <- c(tn, paste0('dcor(', vars[j], ', past)'))
    targets[[length(targets) + 1L]] <- .ctAdequacyCentredDistance(E[use, j, drop = FALSE]^2)
    tn <- c(tn, paste0('dcor(', vars[j], '^2, past)'))
  }
  va <- mean(A * A)
  dcor <- function(B, idx = NULL) {
    Bp <- if (is.null(idx)) B else B[idx, idx]
    sqrt(max(0, mean(A * Bp)) / sqrt(va * mean(B * B)))
  }
  observed <- vapply(targets, dcor, 0)
  m <- length(use)
  perm <- vapply(seq_len(nperm), function(b) {
    idx <- sample.int(m)
    vapply(targets, dcor, 0, idx = idx)
  }, observed)
  perm <- matrix(perm, nrow = length(targets))
  pstat <- (1 + rowSums(perm >= observed)) / (1 + nperm)
  zp <- function(x) stats::qnorm(pmax(x, 1e-300) / 2, lower.tail = FALSE)
  ratio <- observed / rowMeans(perm)
  maxperm <- apply(perm / rowMeans(perm), 2, max)
  w <- which.max(ratio)
  list(statistics = data.frame(family = 'dependence', statistic = tn, observed = observed,
    z = zp(pstat), p = pstat, stringsAsFactors = FALSE),
    family = data.frame(family = 'dependence', statistics = length(targets), worst = tn[w],
      worst_z = zp(pstat[w]), p = (1 + sum(maxperm >= max(ratio))) / (1 + nperm),
      stringsAsFactors = FALSE))
}

# The families compared with generated data. Each statistic is standardised
# by its generated mean and sd and referred to a t on the number of datasets;
# a family's p-value is its statistics' Cauchy combination.
.ctAdequacyGenerated <- function(fit, covariates, which, nsamples, lags, cores) {
  if (is.null(fit$generated) || dim(fit$generated$Y)[1] < nsamples) {
    fit <- suppressMessages(ctGenerateFromFit(fit, nsamples = nsamples, cores = cores))
  }
  model <- .ctFitModelObject(fit)
  vars <- model$manifestNames
  long <- as.data.frame(.ctFitLongData(fit))
  id <- long[[model$subjectIDname]]
  Y <- as.matrix(long[, vars, drop = FALSE])
  G <- fit$generated$Y
  if (dim(G)[2] != nrow(Y)) stop('The generated data do not have the data\'s rows.',
    call. = FALSE)
  G <- G[seq_len(nsamples), , match(vars, dimnames(G)[[3]]), drop = FALSE]
  covs <- if ('covariates' %in% which) .ctScoreCovariateValues(fit, covariates) else NULL
  ids <- unique(id)
  statfun <- function(Ymat) .ctAdequacyStats(Ymat, id, ids, vars, covs)
  observed <- statfun(Y)
  generated <- vapply(seq_len(nsamples), function(b) statfun(matrix(G[b, , ],
    nrow = nrow(Y), dimnames = list(NULL, vars))), observed)
  generated <- matrix(generated, nrow = length(observed), dimnames = list(names(observed), NULL))
  if ('lagged' %in% which) {
    cc <- suppressMessages(ctFitCheckCov(fit, plot = FALSE, keep = 'samples',
      lags = lags, nsamples = nsamples, cores = cores))
    samples <- attr(cc, 'samples')
    cc <- cc[!(cc$lag == 0 & cc$rowvar == cc$colvar), ]
    cc <- cc[!(cc$lag == 0 & cc$rowvar > cc$colvar), ]
    label <- paste0('lagged:', cc$rowvar, '~', cc$colvar, '@', cc$lag)
    slab <- paste0('lagged:', samples$rowvar, '~', samples$colvar, '@', samples$lag)
    lagged <- vapply(seq_len(nsamples), function(b) {
      s <- samples[samples$sample == b, ]
      s$value[match(label, slab[samples$sample == b])]
    }, numeric(length(label)))
    observed <- c(setNames(cc$empirical, label), observed)
    generated <- rbind(matrix(lagged, nrow = length(label), dimnames = list(label, NULL)),
      generated)
  }
  family <- sub(':.*', '', names(observed))
  keep <- family %in% which & is.finite(observed) & apply(generated, 1, function(g)
    sum(is.finite(g)) > 10 && stats::sd(g, na.rm = TRUE) > 0)
  observed <- observed[keep]; generated <- generated[keep, , drop = FALSE]
  family <- family[keep]
  m <- rowMeans(generated, na.rm = TRUE)
  s <- apply(generated, 1, stats::sd, na.rm = TRUE)
  ng <- rowSums(is.finite(generated))
  z <- (observed - m) / s
  pstat <- 2 * stats::pt(-abs(z) / sqrt(1 + 1 / ng), df = ng - 1)
  statistics <- data.frame(family = family, statistic = sub('^[^:]*:', '', names(observed)),
    observed = observed, z = z, p = pstat, stringsAsFactors = FALSE, row.names = NULL)
  families <- do.call(rbind, lapply(unique(family), function(f) {
    rows <- which(family == f)
    w <- rows[which.min(pstat[rows])]
    data.frame(family = f, statistics = length(rows), worst = statistics$statistic[w],
      worst_z = z[w], p = .ctCauchyCombine(pstat[rows]), stringsAsFactors = FALSE)
  }))
  list(fit = fit, statistics = statistics, families = families)
}

# The subject-level, nonlinearity and covariate statistics of one dataset
# (rows in the data's order, `id` per row).
.ctAdequacyStats <- function(Y, id, ids, vars, covs) {
  p <- length(vars)
  rows <- split(seq_len(nrow(Y)), factor(id, levels = ids))
  lag1 <- function(x, y) {
    n <- length(x)
    if (n < 4) return(NA_real_)
    a <- x[-n]; b <- y[-1]
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 3 || stats::sd(a[ok]) == 0 || stats::sd(b[ok]) == 0) return(NA_real_)
    stats::cor(a[ok], b[ok])
  }
  per <- lapply(rows, function(r) {
    y <- Y[r, , drop = FALSE]
    out <- c(colMeans(y, na.rm = TRUE), apply(y, 2, stats::sd, na.rm = TRUE))
    names(out) <- c(paste0('mean_', vars), paste0('sd_', vars))
    for (j in seq_len(p)) for (k in seq_len(p))
      out[paste0('lag1_', vars[j], '_', vars[k])] <- lag1(y[, j], y[, k])
    out
  })
  per <- do.call(rbind, per)
  out <- c()
  for (nm in colnames(per)) out[paste0('subjects:sd(', nm, ')')] <- stats::sd(per[, nm], na.rm = TRUE)
  # Pooled within-subject change against products of deviations.
  dev <- Y; change <- Y; change[] <- NA_real_
  for (r in rows) {
    if (length(r) < 2) next
    dev[r, ] <- sweep(Y[r, , drop = FALSE], 2, colMeans(Y[r, , drop = FALSE], na.rm = TRUE))
    change[r[-length(r)], ] <- Y[r[-1], , drop = FALSE] - Y[r[-length(r)], , drop = FALSE]
  }
  pc <- function(a, b) {
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 10) return(NA_real_)
    stats::cor(a[ok], b[ok])
  }
  for (j in seq_len(p)) {
    for (k in seq_len(p)) for (m in k:p) {
      out[paste0('nonlinear:cor(d', vars[j], ',', vars[k], '*', vars[m], ')')] <-
        pc(change[, j], dev[, k] * dev[, m])
    }
    for (k in seq_len(p)) out[paste0('nonlinear:cor(d', vars[j], '^2,', vars[k], ')')] <-
      pc(change[, j]^2, dev[, k])
  }
  if (length(covs)) for (cn in names(covs)) {
    z <- suppressWarnings(as.numeric(if (is.numeric(covs[[cn]])) covs[[cn]] else
      as.integer(as.factor(covs[[cn]]))))
    for (nm in colnames(per)) out[paste0('covariates:cor(', cn, ',', nm, ')')] <- pc(z, per[, nm])
  }
  out
}

#' @export
print.ctFitAdequacy <- function(x, digits = 3, ...) {
  cat("Adequacy: ", if (x$adequate) "adequate" else "NOT adequate", " at alpha = ",
    x$alpha, if (!x$adequate && isTRUE(x$close_fit)) paste0(
      ", but a close fit (remaining correlations not above ", x$close, ")"), "\n", sep = "")
  notes <- c(if (!is.na(x$projection)) paste0("estimation removed by: ", x$projection),
    if (!is.na(x$nsamples)) paste0(x$nsamples, " generated datasets"))
  if (length(notes)) cat("(", paste(notes, collapse = "; "), ")\n", sep = "")
  cat("\n")
  f <- x$families
  f[c('worst_z', 'p', 'p_adj')] <- lapply(f[c('worst_z', 'p', 'p_adj')], signif, digits)
  print(f, row.names = FALSE)
  if (!is.null(x$size) && nrow(x$size)) {
    cat("\nSize of what remains (root mean square, unbiased):\n")
    s <- x$size[, c('measure', 'estimate', 'lower', 'upper', 'close', 'p_close')]
    s[c('estimate', 'lower', 'upper', 'p_close')] <-
      lapply(s[c('estimate', 'lower', 'upper', 'p_close')], signif, digits)
    print(s, row.names = FALSE)
  }
  invisible(x)
}

# Anderson-Darling against the uniform, for values already put through the
# reference distribution, with Marsaglia and Marsaglia's (2004) limiting
# distribution: c(statistic, p).
.ctAdequacyAD <- function(u) {
  u <- sort(u[is.finite(u)])
  n <- length(u)
  if (n < 10) return(NULL)
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  i <- seq_len(n)
  a2 <- -n - mean((2 * i - 1) * (log(u) + log(1 - rev(u))))
  cdf <- if (a2 < 2) exp(-1.2337141 / a2) / sqrt(a2) * (2.00012 + (0.247105 -
    (0.0649821 - (0.0347962 - (0.011672 - 0.00168691 * a2) * a2) * a2) * a2) * a2) else
    exp(-exp(1.0776 - (2.30695 - (0.43424 - (0.082433 - (0.008056 - 0.0003146 * a2) *
      a2) * a2) * a2) * a2))
  c(a2, min(1, max(0, 1 - cdf)))
}

# The double-centred Euclidean distance matrix of the rows of X, for the
# distance correlation.
.ctAdequacyCentredDistance <- function(X) {
  D <- as.matrix(stats::dist(X))
  D - outer(rowMeans(D), colMeans(D), '+') + mean(D)
}
