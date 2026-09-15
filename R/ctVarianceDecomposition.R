# Variance decomposition -----------------------------------------------------
#
# How much of an indicator's variance is between people, how much is the
# process moving within a person, and how much is measurement error.
#
# The decomposition is a nested law of total variance over the population of
# (person, occasion) pairs the fit was built on. Writing g for the model's
# expected observation given the latent state -- g = E[y | eta] -- and taking
# the expectations over persons i, over each person's own observation times t,
# and over the latent path:
#
#   Var(y) = E[Var(y | eta)]                  measurement
#          + E_i E_t Var_path(g)              process, stochastic
#          + E_i Var_t(E_path[g])             process, deterministic
#          + Var_i(E_{t,path}[g])             between person
#
# Each term is non-negative and the four sum exactly to the total, because each
# split is one application of Var(X) = E[Var(X|Z)] + Var(E[X|Z]). The first
# three are within-person.
#
# Why the within-process term is split in two. A within-person number that
# lumps them together cannot be read: for a stationary process it is all
# diffusion, and for a process with a trend, a time dependent predictor, or one
# still relaxing from its starting distribution, part of it is the mean path
# moving and has nothing to do with the system noise. Separating them is free
# once the moment recursion is running, and it is the only place where
# non-stationarity shows up as a number rather than as a caveat.
#
# What the moments come from. For a linear model the per-row marginal moments
# of the latent state given a person's parameters -- no data, so this is the
# model's own claim rather than a smoothed estimate -- follow from the ordinary
# forward recursion over that person's actual observation times:
#
#   m_1 = T0MEANS,  P_1 = T0VAR
#   m_t = A m_{t-1} + b,  P_t = A P_{t-1} A' + Q
#
# with A, b and Q from the interval. Unequal spacing, missingness, a time
# dependent predictor and a process that has not reached stationarity are then
# all handled exactly rather than assumed away. Collapsing to the stationary
# case -- P -> asymDIFFUSIONcov, m -> asymCINT -- recovers the familiar
# LAMBDA %*% asymDIFFUSIONcov %*% t(LAMBDA) formula, and
# test-ctVarianceDecomposition.R checks that it does as the window grows.
#
# This is exact for a model whose dynamics are linear. It is not exact, and not
# reported, for a model with state-dependent DRIFT, DIFFUSION or LAMBDA cells:
# there the recursion would be a linearisation at a point nothing here chose,
# which is the failure mode CLAUDE.md's "a reported matrix must say where it was
# evaluated" section exists to stop. Such a fit is refused by name, through the
# same detector every other reporting function uses.
#
# Both backends are served from one body. Everything read here --
# `ctExtract(subjectMatrices=TRUE)`, the design accessors in
# R/ctBackendKalman.R, `.ctFitConditionalCells()` -- already means the same
# thing on each.


# Gauss-Hermite nodes and weights for a standard normal ------------------------
#
# Golub-Welsch: the probabilists' Hermite recurrence has zero diagonal and
# off-diagonal sqrt(k), so the nodes are that symmetric tridiagonal matrix's
# eigenvalues and the weights the squared first components of its eigenvectors.
# The measure is already normalised, so the weights sum to one and
# sum(weight * f(node)) is E[f(Z)] with Z standard normal.
#
# Written out rather than taken from statmod, which ctsem does not depend on.
.ctVarDecompGaussHermite <- function(n = 21L) {
  n <- as.integer(n)
  if (n < 3L) stop('quadpoints must be at least 3.', call. = FALSE)
  k <- seq_len(n - 1L)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(k, k + 1L)] <- sqrt(k)
  jacobi[cbind(k + 1L, k)] <- sqrt(k)
  e <- eigen(jacobi, symmetric = TRUE)
  ord <- order(e$values)
  list(node = e$values[ord], weight = (e$vectors[1L, ord])^2)
}


# One interval's transition, intercept and innovation covariance --------------
#
# The intercept uses the block exponential rather than ctsem's own
# `asymCINT + expm(DRIFT dt) (state - asymCINT)` form. The two agree wherever
# both are defined -- expm([[A, c], [0, 0]] dt) has the integral
# int_0^dt expm(A s) ds %*% c in its top right block -- and this one is also
# defined when DRIFT is singular, where asymCINT is not.
#
# The innovation covariance is `.ctNetworkInnovation()` (R/ctGraph.R), Van
# Loan's block exponential, for the same reason: it needs no stationarity
# assumption, where the engine's `asymDIFFUSIONcov - eJAx asymDIFFUSIONcov eJAx'`
# does.
.ctVarDecompStep <- function(drift, diffusioncov, cint, dt, continuoustime) {
  n <- nrow(drift)
  if (continuoustime) {
    block <- rbind(cbind(drift, cint), matrix(0, 1L, n + 1L)) * dt
    e <- as.matrix(expm::expm(block))
    transition <- e[seq_len(n), seq_len(n), drop = FALSE]
    intercept <- e[seq_len(n), n + 1L, drop = FALSE]
  } else {
    steps <- as.integer(round(dt))
    if (steps < 1L) stop(call. = FALSE,
      'A discrete time model needs whole number time steps; found an interval of ', dt, '.')
    transition <- diag(n)
    intercept <- matrix(0, n, 1L)
    for (s in seq_len(steps)) {
      intercept <- drift %*% intercept + cint
      transition <- drift %*% transition
    }
  }
  list(transition = transition, intercept = intercept,
    innovation = .ctNetworkInnovation(drift, diffusioncov, dt, continuoustime))
}


# The marginal moments of one person's trajectory over their own rows ---------
#
# No data enters: these are the model's moments given that person's parameters,
# which is what a decomposition of the model's implied variance needs. A
# smoothed trajectory would answer a different question and would be shrunk
# toward the observations.
#
# Steps are cached by interval. Most designs have a handful of distinct
# intervals and each one costs two matrix exponentials.
.ctVarDecompPersonMoments <- function(mats, times, tdpreds, continuoustime) {
  nlatent <- nrow(mats$DRIFT)
  nmanifest <- nrow(mats$LAMBDA)
  nrows <- length(times)
  latentmean <- matrix(0, nrows, nlatent)
  latentvar <- matrix(0, nrows, nlatent)
  linearmean <- matrix(0, nrows, nmanifest)
  linearvar <- matrix(0, nrows, nmanifest)
  state <- as.numeric(mats$T0MEANS)
  cov <- mats$T0cov
  cache <- list()
  for (ri in seq_len(nrows)) {
    if (ri > 1L) {
      dt <- times[ri] - times[ri - 1L]
      key <- format(dt, digits = 15L)
      step <- cache[[key]]
      if (is.null(step)) {
        step <- .ctVarDecompStep(mats$DRIFT, mats$DIFFUSIONcov, mats$CINT, dt,
          continuoustime)
        cache[[key]] <- step
      }
      state <- as.numeric(step$transition %*% state + step$intercept)
      cov <- step$transition %*% cov %*% t(step$transition) + step$innovation
    }
    # A time dependent predictor is an instantaneous shift at the observed row,
    # applied after the time update -- matching the filter, which does
    # `state[1:nlatent] += (TDPREDEFFECT * tdpreds[rowi])'` at that point.
    if (!is.null(tdpreds) && !is.null(mats$TDPREDEFFECT)) {
      state <- state + as.numeric(mats$TDPREDEFFECT %*% tdpreds[ri, ])
    }
    cov <- (cov + t(cov)) / 2
    latentmean[ri, ] <- state
    latentvar[ri, ] <- diag(cov)
    linearmean[ri, ] <- as.numeric(mats$LAMBDA %*% state + mats$MANIFESTMEANS)
    linearvar[ri, ] <- diag(mats$LAMBDA %*% cov %*% t(mats$LAMBDA))
  }
  list(latentmean = latentmean, latentvar = latentvar,
    linearmean = linearmean, linearvar = linearvar)
}


# The measurement model's contribution, per row and per variable --------------
#
# Given the linear predictor's marginal mean and variance at a row, return the
# three quantities the decomposition needs:
#
#   expected   E_eta[ E[y | eta] ]
#   varmean    Var_eta( E[y | eta] )    -- the process part, on the reported scale
#   condvar    E_eta[ Var(y | eta) ]    -- the measurement part
#
# On the latent scale E[y|eta] is the linear predictor itself, so `expected` and
# `varmean` are what the recursion already produced and `condvar` is the
# measurement variance: MANIFESTcov for a Gaussian indicator, and pi^2/3 for a
# binary or ordinal one, the logistic variance implied by the threshold form
# P(y <= k | eta) = inv_logit(tau_k - eta).
#
# On the response scale the link is integrated over the row's own normal
# marginal by Gauss-Hermite, which is exact to quadrature error.
.ctVarDecompMeasurement <- function(mean, var, type, manifestvar, scale, gh) {
  if (type == 0L) {
    return(list(expected = mean, varmean = var,
      condvar = rep(manifestvar, length(mean))))
  }
  if (identical(scale, 'latent')) {
    return(list(expected = mean, varmean = var,
      condvar = rep(pi^2 / 3, length(mean))))
  }
  sd <- sqrt(pmax(var, 0))
  linear <- matrix(mean, nrow = length(mean), ncol = length(gh$node)) +
    outer(sd, gh$node)
  p <- stats::plogis(linear)
  expected <- as.numeric(p %*% gh$weight)
  list(expected = expected,
    varmean = pmax(as.numeric((p^2) %*% gh$weight) - expected^2, 0),
    condvar = as.numeric((p * (1 - p)) %*% gh$weight))
}


# Design ----------------------------------------------------------------------

# Time dependent predictor values per row, for whichever backend. The julia
# spec keeps them predictor by row, standata row by predictor.
.ctVarDecompTDpredData <- function(fit) {
  model <- .ctFitModelObject(fit)
  if (!length(model$TDpredNames)) return(NULL)
  values <- if (!.ctFitIsJulia(fit) && !is.null(fit$standata$tdpreds)) {
    as.matrix(fit$standata$tdpreds)
  } else t(as.matrix(.ctBackendSpec(fit)$tdpred_data))
  if (ncol(values) != length(model$TDpredNames)) values <- t(values)
  colnames(values) <- model$TDpredNames
  values
}

.ctVarDecompDesign <- function(fit) {
  subject <- .ctFitRowSubject(fit)
  list(subject = subject, time = .ctFitRowTime(fit),
    tdpreds = .ctVarDecompTDpredData(fit),
    nsubjects = max(subject))
}


# Person sources --------------------------------------------------------------
#
# A person is a set of model matrices plus the design rows to evaluate them
# over. Two sources, and on a julia fit they share a body: an individually
# varying parameter is carried as a latent state with no drift and no
# diffusion, so a person *is* a carrier vector, and `ctBackendParMatrices()`
# materialises that person's matrices at it -- applying every transform through
# the same engine code the likelihood uses, rather than reimplementing one
# here. persons='model' draws the carrier vector from the population
# distribution (the carrier block of the augmented T0MEANS and T0cov, shifted
# by that person's time independent predictors); persons='estimated' reads each
# subject's own carrier values off the filter instead.
#
# Each person is given an observed subject's design, so the occasions, spacing,
# missingness and time dependent predictor values stay the ones the fit was
# built on. For persons='model' the donor subject is sampled with replacement.

# Trim the augmented state out of one person's matrices.
#
# The subject and population matrices are already over the real latent
# processes -- except the T0 pair, which keeps the carrier rows because that is
# where the random effects live -- so those are cut here. The par-matrices
# route is asked for `trim = FALSE` so that the carrier block can be read
# whole, and is cut by the same rule .ctBackendTrimAugmented() applies.
.ctVarDecompTrim <- function(mats, nlatent, augmented) {
  if (augmented) {
    for (name in intersect(names(mats), c('DRIFT', 'DIFFUSIONcov'))) {
      mats[[name]] <- mats[[name]][seq_len(nlatent), seq_len(nlatent), drop = FALSE]
    }
    for (name in intersect(names(mats), c('CINT', 'TDPREDEFFECT'))) {
      mats[[name]] <- mats[[name]][seq_len(nlatent), , drop = FALSE]
    }
    if (!is.null(mats$LAMBDA)) mats$LAMBDA <- mats$LAMBDA[, seq_len(nlatent), drop = FALSE]
  }
  mats$T0MEANS <- mats$T0MEANS[seq_len(nlatent), , drop = FALSE]
  mats$T0cov <- mats$T0cov[seq_len(nlatent), seq_len(nlatent), drop = FALSE]
  mats
}

# Subject matrices at the point estimate, for whichever backend.
.ctVarDecompSubjectMatrices <- function(fit) {
  if (.ctFitIsJulia(fit)) {
    point <- fit
    point$estimate$rawposterior <- NULL
    return(ctExtract(point, subjectMatrices = TRUE))
  }
  suppressMessages(stan_constrainsamples(sm = fit$stanmodel,
    standata = fit$standata, samples = matrix(fit$stanfit$rawest, nrow = 1L),
    cores = 1L, savescores = FALSE, savesubjectmatrices = TRUE,
    dokalman = TRUE, onlyfirstrow = FALSE, pcovn = 5))
}

.ctVarDecompNeeded <- c('DRIFT', 'DIFFUSIONcov', 'CINT', 'LAMBDA',
  'MANIFESTMEANS', 'MANIFESTcov', 'T0MEANS', 'T0cov', 'TDPREDEFFECT')

# The carrier state indices of an augmented fit, and nothing when there is no
# augmentation.
.ctVarDecompCarrier <- function(fit, nlatent) {
  augmented <- if (.ctFitIsJulia(fit)) as.integer(.ctBackendSpec(fit)$nlatent_augmented)
    else as.integer(fit$standata$nlatentpop)
  if (!length(augmented) || is.na(augmented) || augmented <= nlatent) return(integer())
  (nlatent + 1L):augmented
}

# One person's matrices, materialised at a carrier vector. julia only -- this is
# the route stan cannot take, and the caller has already said so.
.ctVarDecompPersonAt <- function(fit, tipreds, state, nlatent) {
  drawn <- ctBackendParMatrices(fit, tipreds = tipreds, state = state,
    trim = FALSE)
  mats <- list()
  for (name in .ctVarDecompNeeded) {
    if (!is.null(drawn[[name]])) mats[[name]] <- as.matrix(drawn[[name]])
  }
  .ctVarDecompTrim(mats, nlatent, augmented = TRUE)
}

# Persons for a julia fit: one carrier vector each, drawn or estimated.
.ctVarDecompJuliaPersons <- function(fit, design, nlatent, source, npersons,
  subjects) {
  carrier <- .ctVarDecompCarrier(fit, nlatent)
  tipreds <- if (length(.ctFitModelObject(fit)$TIpredNames))
    .ctFitTIpredData(fit) else NULL
  tipredrow <- function(si) if (!is.null(tipreds)) as.numeric(tipreds[si, ]) else NULL

  donors <- if (identical(source, 'model'))
    subjects[sample.int(length(subjects), npersons, replace = TRUE)] else subjects

  # The subject's own carrier values, for persons='estimated'. These are what
  # the filter has learned about that subject, so they are shrunk toward the
  # population mean -- which is the whole difference between the two sources.
  #
  # Only the carrier block of `subj_T0MEANS` is read. Its dynamic rows are the
  # *smoothed* initial state rather than the model's T0MEANS, so starting a
  # marginal recursion there would begin each person at a data-informed point
  # while still carrying the full prior T0 covariance -- counting the initial
  # spread twice, and reporting it as the mean path moving. The dynamic rows
  # come from the materialised matrices below instead.
  estimated <- if (identical(source, 'estimated') && length(carrier)) {
    .ctVarDecompSubjectMatrices(fit)$subj_T0MEANS
  }

  # One population fetch per donor subject actually used: the augmented T0
  # block depends on the parameters and on that subject's predictors, and on
  # nothing else.
  population <- list()
  for (si in unique(donors)) {
    population[[as.character(si)]] <- ctBackendParMatrices(fit,
      tipreds = tipredrow(si), trim = FALSE)
  }

  # With no carrier states every person with the same donor has the same
  # matrices, so they are materialised once per donor rather than once per
  # person -- a model with no random effects would otherwise pay npersons
  # engine round trips to compute one answer npersons times.
  materialised <- list()
  lapply(seq_along(donors), function(index) {
    si <- donors[index]
    base <- population[[as.character(si)]]
    state <- as.numeric(base$T0MEANS)
    if (length(carrier)) {
      if (identical(source, 'model')) {
        root <- .ctVarDecompCholesky(base$T0cov[carrier, carrier, drop = FALSE])
        state[carrier] <- state[carrier] +
          as.numeric(root %*% stats::rnorm(length(carrier)))
      } else state[carrier] <- as.numeric(estimated[1L, si, carrier, 1L])
    } else {
      key <- as.character(si)
      if (is.null(materialised[[key]])) {
        materialised[[key]] <<- .ctVarDecompPersonAt(fit, tipredrow(si), state,
          nlatent)
      }
      return(list(mats = materialised[[key]],
        rows = which(design$subject == si)))
    }
    list(mats = .ctVarDecompPersonAt(fit, tipredrow(si), state, nlatent),
      rows = which(design$subject == si))
  })
}

# Persons for a stan fit: the subject matrices the filter saved, falling back to
# the population matrix.
#
# The fallback is not a guess. Stan allocates a per-subject array for a matrix
# exactly when that matrix's specification varies by subject (the
# `savesubjectmatrices && (sum(whenmat[..]) || statedep[..])` gate in
# R/ctModelWriter.R), so a matrix with no `subj_` entry is one that is the same
# for everyone and `pop_` is its value. Reading `subj_` alone left LAMBDA and
# the measurement matrices missing for every model that does not vary them.
.ctVarDecompStanPersons <- function(fit, design, nlatent, subjects) {
  extracted <- .ctVarDecompSubjectMatrices(fit)
  lapply(subjects, function(si) {
    mats <- list()
    for (name in .ctVarDecompNeeded) {
      value <- extracted[[paste0('subj_', name)]]
      if (!is.null(value)) {
        mats[[name]] <- array(value[1L, si, , ], dim = dim(value)[3:4])
        next
      }
      value <- extracted[[paste0('pop_', name)]]
      if (is.null(value)) next
      mats[[name]] <- array(value[1L, , ], dim = dim(value)[2:3])
    }
    list(mats = .ctVarDecompTrim(mats, nlatent, augmented = FALSE),
      rows = which(design$subject == si))
  })
}

# A square root of a covariance block that may be singular -- a parameter with
# no random effect has a zero row and column, and chol() refuses those.
.ctVarDecompCholesky <- function(cov) {
  cov <- (cov + t(cov)) / 2
  e <- eigen(cov, symmetric = TRUE)
  values <- pmax(e$values, 0)
  e$vectors %*% diag(sqrt(values), nrow = length(values))
}


# Accumulation ----------------------------------------------------------------
#
# Person by person, then across persons. A person contributes the average over
# its own rows, so persons count equally however many observations each has --
# the estimand is over persons, and weighting by row count would make a
# frequently measured subject a bigger share of the population.
#
# The time term uses the population variance over that person's rows, because
# those rows are the design rather than a sample from anything. The between
# term uses the sample variance over persons, which is an estimate of a
# population variance.
.ctVarDecompPopVar <- function(x) {
  if (length(x) < 2L) return(0)
  mean((x - mean(x))^2)
}


# Refusals --------------------------------------------------------------------

.ctVarDecompCheckModel <- function(fit, scale) {
  cells <- .ctFitConditionalCells(fit)
  if (nrow(cells)) {
    stop('Cells of ', paste(.ctContextReportableMatrices(cells), collapse = ', '),
      ' depend on the latent state or a time dependent predictor, so this ',
      'model has no single DRIFT, DIFFUSION or LAMBDA and the moment ',
      'recursion this function runs would be a linearisation at a point ',
      'nothing here chose. The decomposition is not reported rather than ',
      'reported wrongly. Use ctPhasePortrait() or ctStateDependencePlot() to ',
      'see how the dynamics vary over the state space.', call. = FALSE)
  }
  model <- .ctFitModelObject(fit)
  type <- as.integer(model$manifesttype)
  unsupported <- which(type %in% c(3L, 4L))
  if (length(unsupported)) {
    stop('Count and censored indicators are not supported yet: ',
      paste(model$manifestNames[unsupported], collapse = ', '),
      '. Their conditional variance needs the engine\'s own measurement ',
      'integral rather than the logistic one used here.', call. = FALSE)
  }
  if (identical(scale, 'response') && any(type %in% 2L)) {
    stop("scale='response' is not available for an ordinal indicator: the ",
      'variance of the observed values depends on how the categories are ',
      'coded, which is the user\'s choice rather than the model\'s. Ordinal ',
      "indicator(s): ", paste(model$manifestNames[type %in% 2L], collapse = ', '),
      ". Use scale='latent', which decomposes the linear predictor behind the ",
      'cumulative logit.', call. = FALSE)
  }
  invisible(NULL)
}


#' Decompose model implied variance into between person, process and measurement
#'
#' Splits the variance each indicator of a fitted model implies into a between
#' person part, a within person part from the latent process, and a within
#' person part from measurement error -- with the process part further split
#' into what the system noise contributes and what the mean path moving
#' contributes.
#'
#' @param fit fit object as generated by \code{\link{ctFit}}, from either
#'   backend.
#' @param persons Which population the between person variance refers to.
#'   \code{'model'} draws persons from the fitted population distribution of
#'   the random effects, which is the model's own claim about the population
#'   and is not shrunk. \code{'estimated'} uses the subjects in the data at
#'   their estimated (empirical Bayes) matrices, which describes this sample
#'   but attenuates the between person variance. \code{'auto'}, the default, is
#'   \code{'model'} for a \code{backend='julia'} fit and \code{'estimated'}
#'   for a stan one, saying which it used -- only the julia engine can
#'   materialise the model matrices at a drawn set of random effects.
#' @param scale For a binary or ordinal indicator, whether to decompose the
#'   linear predictor behind the logit (\code{'latent'}, the default, with a
#'   measurement variance of \eqn{\pi^2/3}) or the observed response
#'   (\code{'response'}, integrating the link over each row's own normal
#'   marginal). The two coincide for a Gaussian indicator.
#' @param npersons Number of persons to draw when \code{persons='model'}.
#' @param latents If TRUE, also decompose the latent processes, which have no
#'   measurement component.
#' @param quadpoints Gauss-Hermite nodes used when \code{scale='response'}.
#' @param subjects \code{'all'}, or an integer vector of subjects whose designs
#'   the decomposition is computed over.
#'
#' @details The decomposition is a nested law of total variance over the
#'   population of (person, occasion) pairs at the design the fit was built on.
#'   With \eqn{g = E[y | \eta]} the model's expected observation given the
#'   latent state,
#'
#'   \deqn{Var(y) = E[Var(y|\eta)] + E_i E_t Var(g) + E_i Var_t(E[g]) + Var_i(E[g])}
#'
#'   whose four terms are the four columns reported. They are non-negative and
#'   sum to the total exactly.
#'
#'   \code{within.deterministic} is where non-stationarity appears: it is the
#'   variance of a person's own mean path over their observation times, so it
#'   is zero for a stationary process with no time dependent predictors and
#'   positive when there is a trend, a predictor driven shift, or a process
#'   still relaxing from its starting distribution.
#'
#'   The moments come from the forward recursion over each subject's actual
#'   observation times, with no data entering, so unequal spacing, missingness
#'   and time dependent predictors are handled exactly rather than by assuming
#'   stationarity. The result is exact for a model with linear dynamics. A
#'   model with state dependent \code{DRIFT}, \code{DIFFUSION} or
#'   \code{LAMBDA} cells is refused rather than linearised silently.
#'
#' @return A data frame of class \code{ctVarianceDecomposition}, one row per
#'   variable, with columns \code{variable}, \code{type}, \code{between},
#'   \code{within.deterministic}, \code{within.stochastic},
#'   \code{within.measurement}, \code{within}, \code{total} and the
#'   corresponding proportions. The population the between column refers to,
#'   the scale, and the number of persons are recorded as attributes and shown
#'   by \code{print}.
#'
#' @seealso \code{\link{ctSummaryMatrices}} for the matrices this is built
#'   from, and \code{\link{ctStateDependencePlot}} for a model whose dynamics
#'   vary over the state space.
#'
#' @examples
#' \donttest{
#' ctVarianceDecomposition(ctstantestfit)
#' }
#' @export
ctVarianceDecomposition <- function(fit, persons = c('auto', 'model', 'estimated'),
  scale = c('latent', 'response'), npersons = 200L, latents = TRUE,
  quadpoints = 21L, subjects = 'all') {

  if (!inherits(fit, c('ctStanFit', 'ctJuliaFit'))) {
    stop('fit object is not a ctsem fit!', call. = FALSE)
  }
  persons <- match.arg(persons)
  scale <- match.arg(scale)
  .ctVarDecompCheckModel(fit, scale)

  model <- .ctFitModelObject(fit)
  nlatent <- .ctFitNlatent(fit)
  continuoustime <- isTRUE(model$continuoustime)
  manifestNames <- model$manifestNames
  latentNames <- model$latentNames[seq_len(nlatent)]
  type <- as.integer(model$manifesttype)

  if (identical(persons, 'auto')) {
    persons <- if (.ctFitIsJulia(fit)) 'model' else 'estimated'
    if (identical(persons, 'estimated')) {
      message("persons='estimated' for this stan fit: drawing persons from the ",
        'population distribution needs the model matrices materialised at a ',
        'drawn set of random effects, which only the julia engine can do. The ',
        'between person variance below is therefore the spread of the ',
        'estimated subjects, which shrinkage attenuates.')
    }
  }
  if (identical(persons, 'model') && !.ctFitIsJulia(fit)) {
    stop("persons='model' needs a backend='julia' fit: it materialises the ",
      'model matrices at each drawn set of random effects, which stan computed ',
      'once during sampling and cannot recompute. Use persons=\'estimated\' ',
      'here, reading its between person variance as attenuated by shrinkage, ',
      "or refit with backend='julia'.", call. = FALSE)
  }

  design <- .ctVarDecompDesign(fit)
  wanted <- if (identical(subjects, 'all')) seq_len(design$nsubjects) else
    as.integer(subjects)
  wanted <- wanted[wanted %in% unique(design$subject)]
  if (!length(wanted)) stop('No rows for the requested subjects.', call. = FALSE)

  people <- if (.ctFitIsJulia(fit)) {
    .ctVarDecompJuliaPersons(fit, design, nlatent, persons,
      as.integer(npersons), wanted)
  } else .ctVarDecompStanPersons(fit, design, nlatent, wanted)
  people <- people[vapply(people, function(p) length(p$rows) > 0L, logical(1L))]
  if (!length(people)) stop('No usable persons.', call. = FALSE)

  gh <- if (identical(scale, 'response')) .ctVarDecompGaussHermite(quadpoints) else NULL

  nmanifest <- length(manifestNames)
  nvar <- nmanifest + if (latents) nlatent else 0L
  personmean <- matrix(NA_real_, length(people), nvar)
  persondet <- matrix(NA_real_, length(people), nvar)
  personstoch <- matrix(NA_real_, length(people), nvar)
  personmeas <- matrix(0, length(people), nvar)

  for (pi in seq_along(people)) {
    person <- people[[pi]]
    rows <- person$rows
    moments <- .ctVarDecompPersonMoments(person$mats, design$time[rows],
      if (!is.null(design$tdpreds)) design$tdpreds[rows, , drop = FALSE] else NULL,
      continuoustime)
    for (vi in seq_len(nmanifest)) {
      measured <- .ctVarDecompMeasurement(moments$linearmean[, vi],
        moments$linearvar[, vi], type[vi], person$mats$MANIFESTcov[vi, vi],
        scale, gh)
      personmean[pi, vi] <- mean(measured$expected)
      persondet[pi, vi] <- .ctVarDecompPopVar(measured$expected)
      personstoch[pi, vi] <- mean(measured$varmean)
      personmeas[pi, vi] <- mean(measured$condvar)
    }
    if (latents) for (li in seq_len(nlatent)) {
      vi <- nmanifest + li
      personmean[pi, vi] <- mean(moments$latentmean[, li])
      persondet[pi, vi] <- .ctVarDecompPopVar(moments$latentmean[, li])
      personstoch[pi, vi] <- mean(moments$latentvar[, li])
    }
  }

  out <- data.frame(
    variable = c(manifestNames, if (latents) latentNames),
    type = c(rep('manifest', nmanifest), if (latents) rep('latent', nlatent)),
    between = if (length(people) > 1L) apply(personmean, 2L, stats::var) else rep(0, nvar),
    within.deterministic = colMeans(persondet),
    within.stochastic = colMeans(personstoch),
    within.measurement = colMeans(personmeas),
    stringsAsFactors = FALSE)
  out$within <- out$within.deterministic + out$within.stochastic + out$within.measurement
  out$total <- out$between + out$within
  for (part in c('between', 'within.deterministic', 'within.stochastic',
    'within.measurement', 'within')) {
    out[[paste0('prop.', part)]] <- out[[part]] / out$total
  }
  rownames(out) <- NULL

  attr(out, 'persons') <- persons
  attr(out, 'scale') <- scale
  attr(out, 'npersons') <- length(people)
  attr(out, 'continuoustime') <- continuoustime
  class(out) <- c('ctVarianceDecomposition', 'data.frame')
  out
}


#' @export
print.ctVarianceDecomposition <- function(x, digits = 3L, ...) {
  cat('Model implied variance decomposition\n\n')
  table <- as.data.frame(x)
  show <- c('variable', 'type', 'between', 'within.deterministic',
    'within.stochastic', 'within.measurement', 'total')
  numeric <- vapply(table[show], is.numeric, logical(1L))
  table[show][numeric] <- lapply(table[show][numeric], round, digits)
  print(table[, show], row.names = FALSE)
  cat('\nProportions\n')
  props <- data.frame(variable = x$variable,
    between = round(x$prop.between, digits),
    within = round(x$prop.within, digits),
    of.which.measurement = round(x$within.measurement / x$within, digits))
  print(props, row.names = FALSE)
  cat('\nBetween person variance refers to ',
    if (identical(attr(x, 'persons'), 'model'))
      paste0(attr(x, 'npersons'), ' persons drawn from the fitted population distribution')
    else paste0('the ', attr(x, 'npersons'),
      ' estimated subjects, attenuated by shrinkage'),
    '.\n', sep = '')
  cat('Within person variance is over each person\'s own observation times; ',
    'the deterministic part is the mean path moving.\n', sep = '')
  if (identical(attr(x, 'scale'), 'latent')) {
    cat("Non-Gaussian indicators are decomposed on the latent response scale ",
      '(measurement variance pi^2/3).\n', sep = '')
  } else cat('Non-Gaussian indicators are decomposed on the response scale.\n')
  invisible(x)
}
