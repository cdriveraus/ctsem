# Reporting for a `poprank` fit.
#
# The design constraint is that everything downstream keeps one code path. The
# summary's `popsd` and `rawpopcorr` are already *derived* from the population
# covariance rather than being free parameters -- `.ctBackendAugmentedPopulation()`
# reads `T0cov` and `.ctBackendRandomEffectLevel()` integrates each parameter's
# transform over its own marginal spread by quadrature -- so a reduced-rank fit
# only has to hand them the same thing: one entry per varying parameter, with a
# marginal raw sd and the correlations among them.
#
# So this file computes the full `k x k` population covariance implied by the
# basis block and the regression coefficients, and returns it in exactly the
# shape the full-rank path returns. The single branch is in
# `.ctBackendPopulation()`; nothing else in the summary knows this feature
# exists.

# Raw parameter numbers for a set of labels, from the prepared parameter table.
.ctBackendPopParnumber <- function(spec, labels) {
  table <- as.data.frame(spec$parameter_table, stringsAsFactors = FALSE)
  vapply(labels, function(l) {
    index <- which(!is.na(table$param) & table$param %in% l &
        !is.na(table$parnumber))
    if (!length(index)) NA_integer_ else as.integer(table$parnumber[index[1L]])
  }, integer(1L))
}

# The population covariance of the basis effects, on the raw parameter scale.
#
# `T0cov` is in state units, and `.ctJuliaAugmentRandomEffects()` folds each
# carrier state's `multiplier*meanscale` into the sd transform, so both indices
# have to be divided back out -- the same `scale` the full-rank path divides out
# of its `rawsd`, applied to a covariance rather than a variance.
.ctBackendPopBasisCovariance <- function(t0cov, states, scale, ndraws) {
  r <- length(states)
  out <- array(0, dim = c(ndraws, r, r))
  for (i in seq_len(r)) for (j in seq_len(r)) {
    out[, i, j] <- t0cov[, states[i], states[j]] / (scale[i] * scale[j])
  }
  out
}

# The population moments a `poprank` fit reports: one entry per varying
# parameter, basis effects first, with the regressed effects' variances and
# covariances implied by `Sigma = [[S, S b'], [b S, b S b']]`.
#
# A regressed effect's `parnumber` is its *mean* parameter, which is what the
# quadrature in `.ctBackendRandomEffectLevel()` has to displace to integrate
# that parameter's transform over its population distribution. Its `rawsd` is
# `sqrt((b S b')[i,i])`, which is a derived quantity with a genuine posterior
# spread because `b` and `S` both have one -- so the interval it reports is
# real, unlike the full-rank route's interval on a coordinate the likelihood
# cannot distinguish.
.ctBackendPopRegressionPopulation <- function(spec, samples, layout, flat) {
  regression <- spec$model$popregression
  augmented <- .ctBackendAugmentedSds(spec)
  if (is.null(regression) || is.null(augmented)) return(NULL)
  coefficients <- regression$coefficients
  if (is.null(coefficients) || !nrow(coefficients)) return(NULL)

  sds <- augmented$sds
  basis <- regression$basis
  regressed <- unique(coefficients$param)
  ndraws <- nrow(samples)

  # Basis effects in the order their carrier states appear, which is the order
  # the full-rank path reports and the order `sds` is already in.
  basisparam <- as.character(augmented$param)
  basisorder <- match(basisparam, basis)
  if (anyNA(basisorder)) return(NULL)
  scale <- if (is.null(sds$scale)) rep(1, nrow(sds)) else as.numeric(sds$scale)
  t0cov <- .ctBackendReshape(flat, layout, match("T0cov", layout$matrix))
  S <- .ctBackendPopBasisCovariance(t0cov, sds$row, scale, ndraws)

  # beta, as [regressed, basis] per draw, with columns matched to the basis
  # order S is in.
  betanumber <- matrix(NA_integer_, length(regressed), length(basis),
    dimnames = list(regressed, basis))
  for (ri in seq_len(nrow(coefficients))) {
    betanumber[coefficients$param[ri], coefficients$basis[ri]] <-
      .ctBackendPopParnumber(spec, coefficients$coefficient[ri])
  }
  betanumber <- betanumber[, basisparam, drop = FALSE]
  if (anyNA(betanumber)) return(NULL)
  meannumber <- .ctBackendPopParnumber(spec, regressed)
  if (anyNA(meannumber)) return(NULL)

  k <- length(basis) + length(regressed)
  covariance <- array(0, dim = c(ndraws, k, k))
  ia <- seq_along(basis); ib <- length(basis) + seq_along(regressed)
  covariance[, ia, ia] <- S
  for (d in seq_len(ndraws)) {
    b <- matrix(samples[d, as.integer(betanumber)], length(regressed),
      length(basis))
    Sd <- matrix(S[d, , ], length(basis), length(basis))
    covariance[d, ia, ib] <- Sd %*% t(b)
    covariance[d, ib, ia] <- b %*% Sd
    covariance[d, ib, ib] <- b %*% Sd %*% t(b)
  }

  parnumber <- c(as.integer(augmented$parnumber), as.integer(meannumber))
  parname <- c(.ctBackendParamLabel(augmented$param, augmented$parnumber),
    regressed)
  rawsd <- matrix(vapply(seq_len(k), function(i) sqrt(pmax(covariance[, i, i], 0)),
    numeric(ndraws)), nrow = ndraws)
  rawcorr <- NULL
  if (k > 1L) {
    lower <- which(lower.tri(diag(k)), arr.ind = TRUE)
    rawcorr <- matrix(vapply(seq_len(nrow(lower)), function(entry) {
      i <- lower[entry, 1L]; j <- lower[entry, 2L]
      denominator <- sqrt(covariance[, i, i] * covariance[, j, j])
      ifelse(denominator > 0, covariance[, i, j] / denominator, NA_real_)
    }, numeric(ndraws)), nrow = ndraws)
  }
  list(parnumber = parnumber, param = parname, rawsd = rawsd,
    rawcorr = rawcorr, level = spec$model$subjectIDname,
    covariance = covariance, rank = regression$rank,
    regressed = regressed, approximate = isTRUE(regression$approximate))
}

# The one branch. Everything else in the summary asks for this.
.ctBackendPopulation <- function(spec, samples, layout, flat) {
  out <- .ctBackendPopRegressionPopulation(spec, samples, layout, flat)
  if (!is.null(out)) return(out)
  .ctBackendAugmentedPopulation(spec, samples, layout, flat)
}

# The estimated coefficients, as their own table. This is the object a
# `poprank` fit actually estimated, so it is reported rather than only the
# moments derived from it: `beta[i,j]` says how much of regressed effect `i` a
# one raw unit difference in basis effect `j` implies.
.ctBackendPopRegressionTable <- function(fit, spec, samples, digits = 3) {
  regression <- spec$model$popregression
  if (is.null(regression) || is.null(regression$coefficients)) return(NULL)
  coefficients <- regression$coefficients
  numbers <- .ctBackendPopParnumber(spec, coefficients$coefficient)
  if (anyNA(numbers)) return(NULL)
  draws <- samples[, as.integer(numbers), drop = FALSE]
  out <- data.frame(
    param = coefficients$param, on = coefficients$basis,
    mean = apply(draws, 2L, mean), sd = apply(draws, 2L, stats::sd),
    `2.5%` = apply(draws, 2L, stats::quantile, probs = .025),
    `50%` = apply(draws, 2L, stats::median),
    `97.5%` = apply(draws, 2L, stats::quantile, probs = .975),
    row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
  numeric <- vapply(out, is.numeric, logical(1L))
  out[numeric] <- lapply(out[numeric], round, digits)
  attr(out, 'rank') <- regression$rank
  attr(out, 'approximate') <- isTRUE(regression$approximate)
  attr(out, 'basis') <- regression$basis
  out
}

# What to say under the tables. A reduced-rank fit reports standard deviations
# and correlations it did not estimate, and the reader has to know which.
.ctBackendPopRegressionNote <- function(spec) {
  regression <- spec$model$popregression
  if (is.null(regression)) return(NULL)
  regressed <- paste(unique(regression$coefficients$param), collapse = ', ')
  note <- paste0('poprank = ', regression$rank,
    ': the population sd and correlations of ', regressed,
    ' are implied by their regression on ',
    paste(regression$basis, collapse = ', '),
    ' rather than separately estimated, and each has no variance independent ',
    'of them.')
  if (isTRUE(regression$approximate)) {
    note <- paste0(note, ' This rank is below the ', regression$nmean,
      ' dimensions the model identifies, so it is an approximation: the ',
      'retained parameters absorb what the dropped dimensions carried.')
  }
  note
}
