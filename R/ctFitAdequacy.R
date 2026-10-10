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
# Calibration: each statistic is standardised by its generated mean and sd,
# and a family's p-value is the share of generated datasets whose largest
# |z| over the family's statistics is at least the data's -- the max takes
# care of the many statistics in a family and of their correlation. The
# generated datasets are the fit's own (`ctGenerateFromFit`), at the point
# estimate unless the fit carries draws, which leaves the check slightly
# conservative for what the parameters were fitted to.

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
#'   \item{residuals}{Whether the standardised one-step-ahead prediction
#'     errors are white, as they are under the model whatever the spacing of
#'     the observations: their autocorrelations and cross-correlations within
#'     subject at lags 1 to \code{max(lags)}, and the lag-one autocorrelation
#'     of their squares (volatility that clusters). This is the residual ACF and
#'     CCF, and needs no generated data: its reference is the normal one of a
#'     correlation between independent errors, with a Sidak correction over
#'     the family.}
#'   \item{covariates}{With \code{covariates}: the correlation of each subject
#'     statistic above with each covariate, against the model, whose generated
#'     subjects do not depend on it.}
#' }
#'
#' A family's p-value is the share of generated datasets whose largest
#' standardised discrepancy over the family's statistics is at least the
#' data's, which accounts for the number of statistics and their correlation.
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
  z <- (observed - m) / s
  zgen <- (generated - m) / s
  pstat <- vapply(seq_along(observed), function(k) {
    g <- generated[k, ]; g <- g[is.finite(g)]
    (1 + 2 * min(sum(g <= observed[k]), sum(g >= observed[k]))) / (1 + length(g))
  }, numeric(1))
  family <- sub(':.*', '', names(observed))
  statistics <- data.frame(family = family, statistic = sub('^[^:]*:', '', names(observed)),
    observed = observed, generated_mean = m, generated_sd = s, z = z,
    p = pmin(1, pstat), stringsAsFactors = FALSE, row.names = NULL)
  # The residual family has an analytic reference, so it joins here.
  res <- .ctAdequacyResiduals(fit, max(1L, max(lags)))
  families <- do.call(rbind, lapply(unique(family), function(f) {
    rows <- family == f
    worst <- max(abs(z[rows]))
    ref <- apply(abs(zgen[rows, , drop = FALSE]), 2, max, na.rm = TRUE)
    data.frame(family = f, statistics = sum(rows),
      worst = statistics$statistic[rows][which.max(abs(z[rows]))],
      worst_z = z[rows][which.max(abs(z[rows]))],
      p = (1 + sum(ref >= worst)) / (1 + sum(is.finite(ref))), stringsAsFactors = FALSE)
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

# Whiteness of the standardised prior (one-step-ahead) errors, within subject:
# under the model they are a martingale difference sequence, uncorrelated at
# every lag and in every pair of variables, so a pooled correlation r over n
# pairs is approximately normal with sd 1/sqrt(n) -- the uncentred form, since
# their mean is zero under the model. The squares, less one, test for
# volatility the model does not carry. Sidak over the family.
.ctAdequacyResiduals <- function(fit, maxlag) {
  r <- try(suppressMessages(ctResiduals(fit)), silent = TRUE)
  if (inherits(r, 'try-error')) return(NULL)
  r <- as.data.frame(r)
  vars <- setdiff(names(r), c('Subject', 'Time'))
  r <- r[order(r$Subject, r$Time), , drop = FALSE]
  E <- as.matrix(r[, vars, drop = FALSE])
  subj <- as.integer(factor(r$Subject))
  stat <- list()
  pairs <- function(k) {
    a <- seq_len(nrow(E) - k)
    a[subj[a] == subj[a + k]]
  }
  pc <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    n <- sum(ok)
    if (n < 10) return(c(NA, NA))
    c(sum(x[ok] * y[ok]) / sqrt(sum(x[ok]^2) * sum(y[ok]^2)), n)
  }
  for (k in seq_len(maxlag)) {
    a <- pairs(k)
    for (j in seq_along(vars)) for (l in seq_along(vars)) {
      v <- pc(E[a, j], E[a + k, l])
      stat[[length(stat) + 1L]] <- c(name = paste0(vars[j], '~', vars[l], '@', k), r = v[1], n = v[2])
    }
  }
  a <- pairs(1)
  for (j in seq_along(vars)) {
    v <- pc(E[a, j]^2 - 1, E[a + 1, j]^2 - 1)
    stat[[length(stat) + 1L]] <- c(name = paste0(vars[j], '^2@1'), r = v[1], n = v[2])
  }
  stat <- as.data.frame(do.call(rbind, stat), stringsAsFactors = FALSE)
  stat$r <- as.numeric(stat$r); stat$n <- as.numeric(stat$n)
  stat <- stat[is.finite(stat$r), , drop = FALSE]
  if (!nrow(stat)) return(NULL)
  z <- stat$r * sqrt(stat$n)
  statistics <- data.frame(family = 'residuals', statistic = stat$name,
    observed = stat$r, generated_mean = 0, generated_sd = 1 / sqrt(stat$n), z = z,
    p = 2 * stats::pnorm(-abs(z)), stringsAsFactors = FALSE)
  w <- which.max(abs(z))
  family <- data.frame(family = 'residuals', statistics = nrow(stat),
    worst = stat$name[w], worst_z = z[w],
    p = 1 - (1 - 2 * stats::pnorm(-abs(z[w])))^nrow(stat), stringsAsFactors = FALSE)
  list(statistics = statistics, family = family)
}
