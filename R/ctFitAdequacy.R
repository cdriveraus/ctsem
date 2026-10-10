# Numeric posterior predictive adequacy of a fit.
#
# The graphical checks (ctPostPredPlots, ctFitCheckCov) show what the data and
# the model disagree about; an automated search needs that as numbers with a
# calibrated p-value, and needs checks aimed at the departures it can repair.
# Pooled lagged correlations, the classical ACF/CCF check, cannot see most of
# them: on simulated data a subject-varying drift and a state-dependent drift
# both left every lagged-correlation cell inside its 95 per cent band, while
# the score tests found them at p = 0.0005 and 5e-6 -- averaged over subjects,
# a mixture of regulation rates and a curved drift still decay like one
# exponential. So each family below compares a statistic aimed at one kind of
# departure, computed on the data and on every dataset generated from the fit.
#
# Calibration: each statistic is standardised by its generated mean and sd and
# referred to the normal (as a t on the number of generated datasets), and a
# family takes its smallest p-value with a Sidak correction. An earlier version
# took a family's p-value as the share of generated datasets whose largest |z|
# reached the data's, which handled the correlation between statistics but
# could not go below 1 / (nsamples + 1): with a Holm correction over eight
# families, a covariate effect at z = 5.6 read as adequate. The generated
# datasets are the fit's own (`ctGenerateFromFit`), at the point estimate
# unless the fit carries draws, which leaves the check slightly conservative
# for what the parameters were fitted to.

#' Posterior predictive adequacy of a fitted model, as tests
#'
#' Compares statistics of the data with the same statistics of datasets
#' generated from the fit, in families aimed at the ways a continuous time
#' model can fall short, and gives each family a p-value.
#'
#' \describe{
#'   \item{lagged}{The lagged correlations within and across variables (the
#'     ACF and CCF, lags 0 to \code{max(lags)}), as \code{\link{ctFitCheckCov}}
#'     computes them: wrong dynamics, wrong measurement error, a missing
#'     process.}
#'   \item{subjects}{How much subjects differ in their own mean, standard
#'     deviation and lag-one autocorrelation and cross-correlation of each
#'     variable: individual differences the model does not carry.}
#'   \item{nonlinear}{Whether change depends on level nonlinearly: the
#'     correlation, pooled within subjects, of each variable's change with the
#'     products of the variables' deviations from their subject means, and of
#'     the squared change with those deviations: a drift or a diffusion that
#'     depends on the state.}
#'   \item{covariates}{With \code{covariates}: the correlation of each subject
#'     statistic above with each covariate, against the model, whose generated
#'     subjects do not depend on it.}
#' }
#'
#' Four more families are built on the one-step-ahead prediction errors, each
#' variable's error standardised given the row's other errors, so that it
#' carries the variable's name. Under the model those errors are standard
#' normal and independent of everything before their row, whatever the
#' spacing of the observations, so these families need no generated data:
#' their references are exact for a linear Gaussian model with known
#' parameters (Sidak over each family, or a permutation where stated).
#' \describe{
#'   \item{residuals}{Their autocorrelations and cross-correlations within
#'     subject at lags 1 to \code{max(lags)}, and the lag-one autocorrelation
#'     of their squares (volatility that clusters): the residual ACF and CCF.}
#'   \item{distribution}{Each variable's errors against the standard normal,
#'     and each row's squared whitened distance against the chi-square, by
#'     Anderson-Darling: skew, heavy tails, outlying occasions.}
#'   \item{fit}{Each subject's summed squared distance against the
#'     chi-square on its number of observations: whether subjects are fitted
#'     more unevenly than the model allows, and the worst one.}
#'   \item{dependence}{The distance correlation of each variable's error, and
#'     of its square, with what came before -- the predicted latent states, the
#'     previous errors, the time interval -- by permutation: a conditional mean
#'     or variance the model gets wrong, in any functional form.}
#' }
#'
#' For the families compared with generated data, each statistic is
#' standardised by its mean and standard deviation over the generated
#' datasets and referred to the normal, and a family's p-value is its
#' smallest, Sidak-corrected for the number of statistics -- conservative
#' when they are correlated, as they are. The references ignore that the parameters were estimated from the same
#' data, which makes the checks conservative for what the parameters fit.
#' Detection is not localization: a failure in one variable's equation
#' reaches the others' errors through the fitted covariance.
#' \code{adequate} is \code{TRUE} when no family's Holm-adjusted p-value is
#' below \code{alpha}. The statistics are descriptive and can miss what they
#' do not measure; \code{\link{ctScoreTest}} is the more sensitive check for
#' the departures it covers, and \code{\link{ctModelImprove}} uses both.
#'
#' @param fit A fit from \code{\link{ctFit}}, either backend.
#' @param covariates Optional subject-level covariates, as for
#'   \code{\link{ctScoreTest}}.
#' @param nsamples Number of generated datasets. Reused from \code{fit} when it
#'   already carries at least as many.
#' @param lags Observation lags for the lagged family.
#' @param alpha Level for \code{adequate}.
#' @param cores Cores for generating.
#' @return An object of class \code{ctFitAdequacy}: \code{families} (one row
#'   per family: its p-value, the Holm-adjusted one, the number of
#'   statistics and the worst), \code{statistics} (every statistic: observed,
#'   generated mean and sd, z, two-sided p), \code{adequate}, and \code{fit}
#'   (with its generated data, for reuse).
#' @seealso \code{\link{ctPostPredPlots}} to see the same comparisons.
#' @export
ctFitAdequacy <- function(fit, covariates = NULL, nsamples = 100, lags = 0:3,
  alpha = 0.05, cores = 2) {
  if (!.ctFitIsJulia(fit) && !inherits(fit, 'ctStanFit')) stop('Not a ctsem fit.',
    call. = FALSE)
  if (exists('.Random.seed', envir = globalenv(), inherits = FALSE)) {
    oldseed <- get('.Random.seed', envir = globalenv(), inherits = FALSE)
    on.exit(assign('.Random.seed', oldseed, envir = globalenv()), add = TRUE)
  } else on.exit(rm('.Random.seed', envir = globalenv()), add = TRUE)
  set.seed(20261010L)
  if (is.null(fit$generated) || dim(fit$generated$Y)[1] < nsamples) {
    fit <- suppressMessages(ctGenerateFromFit(fit, nsamples = nsamples, cores = cores))
  }
  model <- .ctFitModelObject(fit)
  vars <- model$manifestNames
  long <- as.data.frame(.ctFitLongData(fit))
  idname <- model$subjectIDname
  id <- long[[idname]]
  Y <- as.matrix(long[, vars, drop = FALSE])
  G <- fit$generated$Y
  if (dim(G)[2] != nrow(Y)) stop('The generated data do not have the data\'s rows.',
    call. = FALSE)
  G <- G[seq_len(nsamples), , match(vars, dimnames(G)[[3]]), drop = FALSE]

  covs <- if (is.null(covariates)) NULL else .ctScoreCovariateValues(fit, covariates)
  ids <- unique(id)
  statfun <- function(Ymat) .ctAdequacyStats(Ymat, id, ids, vars, covs)
  observed <- statfun(Y)
  generated <- vapply(seq_len(nsamples), function(b) statfun(matrix(G[b, , ],
    nrow = nrow(Y), dimnames = list(NULL, vars))), observed)

  # The lagged family, as ctFitCheckCov defines it, from the same datasets.
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
  generated <- rbind(matrix(lagged, nrow = length(label)), generated)
  rownames(generated) <- names(observed)

  keep <- is.finite(observed) & apply(generated, 1, function(g) sum(is.finite(g)) > 10 &&
    stats::sd(g, na.rm = TRUE) > 0)
  observed <- observed[keep]; generated <- generated[keep, , drop = FALSE]
  m <- rowMeans(generated, na.rm = TRUE)
  s <- apply(generated, 1, stats::sd, na.rm = TRUE)
  ng <- rowSums(is.finite(generated))
  z <- (observed - m) / s
  # A new draw against a normal whose mean and sd come from ng draws: the
  # prediction-interval t, so a p-value is not floored at 1 / (nsamples + 1).
  pstat <- 2 * stats::pt(-abs(z) / sqrt(1 + 1 / ng), df = ng - 1)
  family <- sub(':.*', '', names(observed))
  statistics <- data.frame(family = family, statistic = sub('^[^:]*:', '', names(observed)),
    observed = observed, generated_mean = m, generated_sd = s, z = z,
    p = pstat, stringsAsFactors = FALSE, row.names = NULL)
  # The innovation families have exact references, so they join here.
  res <- .ctAdequacyInnovations(fit, max(1L, max(lags)))
  families <- do.call(rbind, lapply(unique(family), function(f) {
    rows <- which(family == f)
    w <- rows[which.min(pstat[rows])]
    data.frame(family = f, statistics = length(rows), worst = statistics$statistic[w],
      worst_z = z[w], p = 1 - (1 - pstat[w])^length(rows), stringsAsFactors = FALSE)
  }))
  if (!is.null(res)) {
    statistics <- rbind(statistics, res$statistics)
    families <- rbind(families, res$family)
  }
  families$p_adj <- stats::p.adjust(families$p, 'holm')
  rownames(families) <- NULL
  structure(list(families = families, statistics = statistics,
    adequate = all(families$p_adj >= alpha), alpha = alpha, nsamples = nsamples,
    fit = fit), class = 'ctFitAdequacy')
}

# The subject-level, nonlinearity, fit-free and covariate statistics of one
# dataset (rows in the data's order, `id` per row).
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
  cat("Posterior predictive adequacy (", x$nsamples, " generated datasets): ",
    if (x$adequate) "adequate" else "NOT adequate", " at alpha = ", x$alpha, "\n\n", sep = "")
  f <- x$families
  f[c('worst_z', 'p', 'p_adj')] <- lapply(f[c('worst_z', 'p', 'p_adj')], signif, digits)
  print(f, row.names = FALSE)
  invisible(x)
}

# The families built on the one-step-ahead prediction errors. Under the model
# the errors of a row, whitened by their predicted covariance, are standard
# normal and independent of everything before the row, whatever the spacing
# of the observations: the prediction-error decomposition is the time-series
# form of the whitening Driver's SEM misfit work (papers/2025 SEM misfit) uses
# on cross-sectional residuals. So these families have exact references and
# need no generated data:
#
#   residuals    lagged correlations of the errors within subject (the
#                residual ACF/CCF), and of their squares at lag one.
#   distribution each variable's conditional error against N(0, 1), and each
#                row's squared whitened distance against chi-square
#                (Anderson-Darling).
#   fit          each subject's summed squared distance against chi-square
#                on its number of observations: the spread of how well the
#                model fits subjects, and the worst subject.
#   dependence   the distance correlation of each variable's conditional
#                error, and of its square, with what came before -- the
#                predicted latent states, the previous row's errors, the time
#                interval: a conditional mean or variance the model gets
#                wrong, in any functional form.
#
# Errors carry variable labels by conditioning rather than by a rotation: each
# variable's error given the row's other errors, standardised (Driver's
# variable-labelled conditional residuals), so a result names a measurement
# equation. Detection is still not localization -- misfit in one equation
# reaches the others through the fitted covariance.
#
# The references ignore that the parameters were estimated from these data,
# which makes each family conservative for what the parameters fit (the
# drift absorbs a linear lag-one dependence, say). The SEM work found a refit
# bootstrap necessary for calibrated levels; here the score tests are the
# estimation-aware layer, and these families are its safety net.
.ctAdequacyInnovations <- function(fit, maxlag, nperm = 299L, nmax = 800L) {
  k <- try(suppressMessages(ctKalmanArray(fit)), silent = TRUE)
  if (inherits(k, 'try-error')) return(NULL)
  vars <- colnames(k$y); if (is.null(vars)) vars <- paste0('Y', seq_len(ncol(k$y)))
  p <- length(vars)
  n <- nrow(k$y)
  id <- as.integer(factor(k$id))
  time <- as.numeric(k$time)
  E <- matrix(NA_real_, n, p, dimnames = list(NULL, vars))
  Q <- rep(NA_real_, n); d <- integer(n)
  for (r in seq_len(n)) {
    obs <- which(!is.na(k$y[r, ]))
    if (!length(obs)) next
    v <- k$errprior[1, r, obs]
    S <- matrix(k$ypriorcov[1, r, obs, obs], length(obs))
    P <- try(solve(S), silent = TRUE)
    if (inherits(P, 'try-error')) next
    Pv <- as.numeric(P %*% v)
    E[r, obs] <- Pv / sqrt(diag(P))
    Q[r] <- sum(v * Pv); d[r] <- length(obs)
  }
  stats <- list(); fams <- list()
  add <- function(family, name, observed, z, p) stats[[length(stats) + 1L]] <<-
    data.frame(family = family, statistic = name, observed = observed,
      generated_mean = NA_real_, generated_sd = NA_real_, z = z, p = p,
      stringsAsFactors = FALSE)
  zp <- function(p) stats::qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)
  family <- function(name, rows, p) {
    w <- which.min(p)
    fams[[length(fams) + 1L]] <<- data.frame(family = name, statistics = length(p),
      worst = rows[w], worst_z = zp(p[w]), p = 1 - (1 - p[w])^length(p),
      stringsAsFactors = FALSE)
  }

  # residuals: lagged correlations within subject.
  pairs <- function(lag) { a <- seq_len(n - lag); a[id[a] == id[a + lag]] }
  # Uncentred for the errors, whose mean is zero under the model; centred for
  # their squares, whose mean is a matter of the error variance, which the
  # distribution family judges.
  pc <- function(x, y, centre = FALSE) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) return(c(NA, NA))
    x <- x[ok]; y <- y[ok]
    if (centre) { x <- x - mean(x); y <- y - mean(y) }
    c(sum(x * y) / sqrt(sum(x^2) * sum(y^2)), sum(ok))
  }
  nm <- character(); pv <- numeric()
  for (lag in seq_len(maxlag)) {
    a <- pairs(lag)
    for (j in seq_len(p)) for (l in seq_len(p)) {
      v <- pc(E[a, j], E[a + lag, l])
      if (!is.finite(v[1])) next
      z <- v[1] * sqrt(v[2]); nm <- c(nm, paste0(vars[j], '~', vars[l], '@', lag))
      pv <- c(pv, 2 * stats::pnorm(-abs(z))); add('residuals', utils::tail(nm, 1), v[1], z, utils::tail(pv, 1))
    }
  }
  a <- pairs(1)
  for (j in seq_len(p)) {
    v <- pc(E[a, j]^2, E[a + 1, j]^2, centre = TRUE)
    if (!is.finite(v[1])) next
    z <- v[1] * sqrt(v[2]); nm <- c(nm, paste0(vars[j], '^2@1'))
    pv <- c(pv, 2 * stats::pnorm(-abs(z))); add('residuals', utils::tail(nm, 1), v[1], z, utils::tail(pv, 1))
  }
  if (length(pv)) family('residuals', nm, pv)

  # distribution: Anderson-Darling against the model's reference.
  nm <- character(); pv <- numeric()
  for (j in seq_len(p)) {
    ad <- .ctAdequacyAD(stats::pnorm(E[, j]))
    if (is.null(ad)) next
    nm <- c(nm, paste0('normal(', vars[j], ')')); pv <- c(pv, ad[2])
    add('distribution', utils::tail(nm, 1), ad[1], zp(ad[2]), ad[2])
  }
  ad <- .ctAdequacyAD(stats::pchisq(Q[d > 0], d[d > 0]))
  if (!is.null(ad)) {
    nm <- c(nm, 'chisq(row distance)'); pv <- c(pv, ad[2])
    add('distribution', 'chisq(row distance)', ad[1], zp(ad[2]), ad[2])
  }
  if (length(pv)) family('distribution', nm, pv)

  # fit: each subject's summed distance.
  Qi <- tapply(Q, id, sum, na.rm = TRUE); di <- tapply(d, id, sum)
  ok <- di > 0
  u <- stats::pchisq(Qi[ok], di[ok])
  ad <- .ctAdequacyAD(u)
  ext <- 2 * pmin(u, 1 - u)
  worst <- which.min(ext)
  pext <- 1 - (1 - ext[worst])^length(u)
  labels <- levels(factor(k$id))[as.integer(names(u))]
  if (!is.null(ad)) {
    add('fit', 'chisq(subject distances)', ad[1], zp(ad[2]), ad[2])
    add('fit', paste0('worst subject: ', labels[worst]), Qi[ok][worst] / di[ok][worst],
      zp(pext), pext)
    family('fit', c('chisq(subject distances)', paste0('worst subject: ', labels[worst])),
      c(ad[2], pext))
  }

  # dependence: distance correlation with the past.
  prev <- c(NA_integer_, seq_len(n - 1L))
  prev[c(TRUE, id[-1] != id[-n])] <- NA_integer_
  eta <- matrix(k$etaprior[1, , ], nrow = n)
  X <- cbind(eta, E[prev, , drop = FALSE], dt = time - time[pmax(prev, 1L)])
  use <- which(!is.na(prev) & stats::complete.cases(X) & stats::complete.cases(E))
  if (length(use) >= 50L) {
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
    for (t in seq_along(targets)) add('dependence', tn[t], observed[t], zp(pstat[t]), pstat[t])
    maxperm <- apply(perm / rowMeans(perm), 2, max)
    w <- which.max(observed / rowMeans(perm))
    fams[[length(fams) + 1L]] <- data.frame(family = 'dependence', statistics = length(targets),
      worst = tn[w], worst_z = zp(pstat[w]),
      p = (1 + sum(maxperm >= max(observed / rowMeans(perm)))) / (1 + nperm),
      stringsAsFactors = FALSE)
  }
  list(statistics = do.call(rbind, stats), family = do.call(rbind, fams))
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
