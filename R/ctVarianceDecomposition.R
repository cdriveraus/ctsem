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
# Messages suppressed, and said once by the print method instead. This is
# called once per person -- up to npersons times -- and each call otherwise
# announces which cells are state dependent and where they were evaluated, so a
# single decomposition printed the same paragraph four hundred times. What that
# paragraph says is still owed to the reader, and the attributes carry it.
.ctVarDecompPersonAt <- function(fit, tipreds, state, nlatent) {
  drawn <- suppressMessages(ctBackendParMatrices(fit, tipreds = tipreds,
    state = state, trim = FALSE))
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
    population[[as.character(si)]] <- suppressMessages(ctBackendParMatrices(fit,
      tipreds = tipredrow(si), trim = FALSE))
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

# The matrices that turn a latent state into an expected observation. A cell of
# one of these that depends on the state is the one kind of nonlinearity the
# simulation route cannot absorb: it draws the *path* from the engine, and then
# needs the measurement model in R to turn each drawn state into an expected
# observation.
.ctVarDecompMeasurementMatrices <- c('LAMBDA', 'MANIFESTMEANS', 'MANIFESTVAR',
  'MANIFESTcov', 'Jy', 'THRESHOLDS')

.ctVarDecompStateDependent <- function(fit, measurementonly = FALSE) {
  cells <- .ctFitConditionalCells(fit)
  if (measurementonly) {
    cells <- cells[cells$matrix %in% .ctVarDecompMeasurementMatrices, , drop = FALSE]
  }
  cells
}

.ctVarDecompCheckModel <- function(fit, method, scale) {
  cells <- .ctVarDecompStateDependent(fit, measurementonly = !identical(method, 'moment'))
  if (nrow(cells)) {
    stop('Cells of ', paste(.ctContextReportableMatrices(cells), collapse = ', '),
      ' depend on the latent state or a time dependent predictor, so they have ',
      'no single value and ',
      if (identical(method, 'moment'))
        paste0('the moment recursion would be a linearisation at a point ',
          "nothing here chose. Use method='simulation', which draws the ",
          'trajectories from the engine instead.')
      else paste0('the simulation route cannot turn a drawn state into an ',
        'expected observation without re-materialising the measurement model ',
        'at every row. Neither route reports a number for this model.'),
      ' ctPhasePortrait() and ctStateDependencePlot() show how the model varies ',
      'over the state space.', call. = FALSE)
  }
  model <- .ctFitModelObject(fit)
  type <- as.integer(model$manifesttype)
  unsupported <- which(type %in% c(3L, 4L))
  if (length(unsupported)) {
    stop('Count and censored indicators are not supported yet: ',
      paste(model$manifestNames[unsupported], collapse = ', '),
      ". Their conditional variance needs the engine's own measurement ",
      'integral rather than the logistic one used here.', call. = FALSE)
  }
  if (identical(scale, 'response') && any(type %in% 2L)) {
    stop("scale='response' is not available for an ordinal indicator: the ",
      'variance of the observed values depends on how the categories are ',
      "coded, which is the user's choice rather than the model's. Ordinal ",
      'indicator(s): ', paste(model$manifestNames[type %in% 2L], collapse = ', '),
      ". Use scale='latent', which decomposes the linear predictor behind the ",
      'cumulative logit.', call. = FALSE)
  }
  invisible(NULL)
}


# The moment route ------------------------------------------------------------
#
# Person by person, then across persons. A person contributes the average over
# its own rows, so persons count equally however many observations each has --
# the estimand is over persons, and weighting by row count would make a
# frequently measured subject a bigger share of the population.
.ctVarDecompMomentComponents <- function(fit, design, people, nlatent, latents,
  type, scale, gh, continuoustime, nmanifest) {
  nvar <- nmanifest + if (latents) nlatent else 0L
  personmean <- matrix(NA_real_, length(people), nvar)
  persondet <- matrix(NA_real_, length(people), nvar)
  personstoch <- matrix(NA_real_, length(people), nvar)
  personmeas <- matrix(0, length(people), nvar)

  for (index in seq_along(people)) {
    person <- people[[index]]
    rows <- person$rows
    moments <- .ctVarDecompPersonMoments(person$mats, design$time[rows],
      if (!is.null(design$tdpreds)) design$tdpreds[rows, , drop = FALSE] else NULL,
      continuoustime)
    for (vi in seq_len(nmanifest)) {
      measured <- .ctVarDecompMeasurement(moments$linearmean[, vi],
        moments$linearvar[, vi], type[vi], person$mats$MANIFESTcov[vi, vi],
        scale, gh)
      personmean[index, vi] <- mean(measured$expected)
      persondet[index, vi] <- .ctVarDecompPopVar(measured$expected)
      personstoch[index, vi] <- mean(measured$varmean)
      personmeas[index, vi] <- mean(measured$condvar)
    }
    if (latents) for (li in seq_len(nlatent)) {
      vi <- nmanifest + li
      personmean[index, vi] <- mean(moments$latentmean[, li])
      persondet[index, vi] <- .ctVarDecompPopVar(moments$latentmean[, li])
      personstoch[index, vi] <- mean(moments$latentvar[, li])
    }
  }

  list(
    between = if (length(people) > 1L) apply(personmean, 2L, stats::var) else rep(0, nvar),
    within.deterministic = colMeans(persondet),
    within.stochastic = colMeans(personstoch),
    within.measurement = colMeans(personmeas),
    npersons = length(people))
}


# Laplace fits: one between column per level ----------------------------------
#
# A Laplace fit's random effects are not carrier states, so the carrier route
# above sees nothing and would report a between person variance of zero for
# every such model -- including every multilevel one, since a level above the
# subject is the thing only Laplace can integrate out.
#
# What it has instead is a level-restricted trajectory. `randomEffects=` names
# a level and means "this level and every level outside it", so for
# `id = c('subject','study')` the three settings give a subject's own matrices,
# the matrices its study shares, and the population's. Each subject therefore
# has a time-averaged expected value at each level, and the differences between
# successive levels are that level's own effect:
#
#     g_subject - g_study   the subject's departure from its study
#     g_study   - g_pop     the study's departure from the population
#
# Those are orthogonal by construction -- each is a difference of successive
# conditional means -- so their variances add, and `between` is their sum with
# one column per level beside it. The variance of a level's deviation is taken
# over the *units of that level*: over subjects for the innermost, over studies
# for the next, which is what makes a study with more subjects in it count once.
#
# These are the fitted modes, so every level is shrunk toward the one outside
# it, and the outer levels hardest because they have the fewest units. That is
# the same caveat persons='estimated' carries on the augmented route and it is
# the only thing on offer here: drawing a person would need a draw from each
# level's own effect covariance, which the engine does not expose.
.ctVarDecompLaplaceLevels <- function(fit) {
  levels <- .ctBackendSpec(fit)$laplace$levels
  vapply(levels, function(x) as.character(x$name)[1L], character(1L))
}

# Subject matrices built from `level` and every level outside it.
.ctVarDecompLevelMatrices <- function(fit, level) {
  scores <- .ctBackendKalmanRaw(fit, fit$estimate$raw, subjectmatrices = TRUE,
    fields = "subject_loglik", randomEffects = level)
  flat <- array(scores$subject_matrices,
    dim = c(1L, dim(scores$subject_matrices)))
  .ctBackendSubjectMatrices(fit, flat)
}

# Which unit of each level every subject belongs to.
#
# From `spec$data` rather than `.ctFitLongData()`, which carries the subject id
# and not the levels above it.
.ctVarDecompLevelUnits <- function(fit, levels) {
  data <- as.data.frame(.ctBackendSpec(fit)$data, stringsAsFactors = FALSE)
  subject <- .ctFitRowSubject(fit)
  first <- match(seq_len(max(subject)), subject)
  out <- lapply(levels, function(name) {
    if (!name %in% names(data)) return(NULL)
    as.character(data[[name]])[first]
  })
  names(out) <- levels
  out
}

# One subject's matrices, pulled out of a level's subject-matrix arrays.
.ctVarDecompLevelPerson <- function(matrices, si, nlatent) {
  mats <- list()
  for (name in .ctVarDecompNeeded) {
    value <- matrices[[paste0('subj_', name)]]
    if (is.null(value)) next
    mats[[name]] <- array(value[1L, si, , ], dim = dim(value)[3:4])
  }
  .ctVarDecompTrim(mats, nlatent, augmented = FALSE)
}

.ctVarDecompLaplaceComponents <- function(fit, design, nlatent, latents, type,
  scale, gh, continuoustime, nmanifest, subjects) {

  levels <- .ctVarDecompLaplaceLevels(fit)
  units <- .ctVarDecompLevelUnits(fit, levels)
  # Innermost first, then each level outside it, then the population -- the
  # order `randomEffects=` means, and the order the differences are taken in.
  steps <- c(levels, 'population')
  matrices <- lapply(steps, function(level) .ctVarDecompLevelMatrices(fit, level))
  names(matrices) <- steps

  nvar <- nmanifest + if (latents) nlatent else 0L
  value <- lapply(steps, function(...) matrix(NA_real_, length(subjects), nvar))
  names(value) <- steps
  persondet <- matrix(NA_real_, length(subjects), nvar)
  personstoch <- matrix(NA_real_, length(subjects), nvar)
  personmeas <- matrix(0, length(subjects), nvar)

  for (index in seq_along(subjects)) {
    si <- subjects[index]
    rows <- which(design$subject == si)
    tdpreds <- if (!is.null(design$tdpreds))
      design$tdpreds[rows, , drop = FALSE] else NULL
    for (level in steps) {
      mats <- .ctVarDecompLevelPerson(matrices[[level]], si, nlatent)
      moments <- .ctVarDecompPersonMoments(mats, design$time[rows], tdpreds,
        continuoustime)
      for (vi in seq_len(nmanifest)) {
        measured <- .ctVarDecompMeasurement(moments$linearmean[, vi],
          moments$linearvar[, vi], type[vi], mats$MANIFESTcov[vi, vi], scale, gh)
        value[[level]][index, vi] <- mean(measured$expected)
        # The within person terms belong to the subject's own trajectory, so
        # they are read from the innermost level and the outer passes only
        # supply the level means.
        if (identical(level, steps[1L])) {
          persondet[index, vi] <- .ctVarDecompPopVar(measured$expected)
          personstoch[index, vi] <- mean(measured$varmean)
          personmeas[index, vi] <- mean(measured$condvar)
        }
      }
      if (latents) for (li in seq_len(nlatent)) {
        vi <- nmanifest + li
        value[[level]][index, vi] <- mean(moments$latentmean[, li])
        if (identical(level, steps[1L])) {
          persondet[index, vi] <- .ctVarDecompPopVar(moments$latentmean[, li])
          personstoch[index, vi] <- mean(moments$latentvar[, li])
        }
      }
    }
  }

  # Each level's own deviation, and its variance over that level's units.
  contribution <- matrix(0, length(levels), nvar,
    dimnames = list(levels, NULL))
  for (l in seq_along(levels)) {
    deviation <- value[[steps[l]]] - value[[steps[l + 1L]]]
    unit <- units[[levels[l]]]
    for (vi in seq_len(nvar)) {
      d <- deviation[, vi]
      if (!is.null(unit)) d <- as.numeric(tapply(d, unit[subjects], mean))
      contribution[l, vi] <- if (length(d) > 1L) stats::var(d) else 0
    }
  }

  list(
    between = colSums(contribution),
    levels = contribution,
    within.deterministic = colMeans(persondet),
    within.stochastic = colMeans(personstoch),
    within.measurement = colMeans(personmeas),
    npersons = length(subjects))
}


# The simulation route --------------------------------------------------------
#
# For a model whose dynamics depend on the state there is no moment recursion to
# run: the transition is a different transition at every point the trajectory
# visits. The engine draws the trajectory instead, through the same state pass
# `ctGenerate(intoverstates = FALSE)` uses, and the decomposition is taken over
# the draws.
#
# The one thing that needs care is that each person must be given *several*
# paths, not one. With a single path per person the average over that person's
# rows is not their mean -- it wanders with the path -- and the wander lands in
# the between person term, inflating it by an amount that depends on how
# autocorrelated the process is. So a person is drawn once and its trajectory
# redrawn `npaths` times:
#
#   `ctsem_state_layout()` says where each subject's innovations start, and the
#   first `nlatent` entries of that block are its initial state draw. An
#   individually varying parameter is a carrier state with no drift and no
#   diffusion, so pinning the carrier entries of that initial block and
#   redrawing everything after them is the same person on a new path.
#
# Pinning the carrier *entries* pins the carrier *values* only when T0VAR has no
# covariance between the dynamic states and the carriers -- the factor applied
# to the draw is triangular, so a cross block would let a redrawn dynamic entry
# move the carrier. That is checked rather than assumed.
#
# Two Monte Carlo corrections, both exact in expectation and both computed from
# the draws themselves. The average over paths estimates a person's mean path
# with error, so
#
#   * the variance over time of that average carries an extra
#     mean_t Var_paths(g_t - gbar) / npaths, and
#   * the variance over persons of a person's overall average carries an extra
#     mean_persons Var_paths(gbar) / npaths.
#
# Without them the deterministic and between terms both grow as npaths falls,
# which is the same finite-path bias in a different place.
.ctVarDecompSimulation <- function(fit, design, nlatent, latents, type, scale,
  npersons, npaths, nmanifest, subjects) {

  if (npaths < 2L) stop('npaths must be at least 2: the simulation route ',
    "separates a person's mean path from the variation around it, and one ",
    'path cannot.', call. = FALSE)
  layout <- .ctBackendStateLayout(fit)
  augmented <- as.integer(layout$nlatent)
  carrier <- if (augmented > nlatent) (nlatent + 1L):augmented else integer()
  population <- suppressMessages(ctBackendParMatrices(fit, trim = FALSE))

  tipreds <- if (length(.ctFitModelObject(fit)$TIpredNames))
    .ctFitTIpredData(fit) else NULL
  raw <- fit$estimate$raw
  nrows <- length(design$subject)
  base <- matrix(0, nmanifest, nrows)

  if (length(carrier)) {
    cross <- population$T0cov[carrier, seq_len(nlatent), drop = FALSE]
    if (any(abs(cross) > 1e-10)) {
      stop('T0VAR has covariance between the latent processes and the ',
        'individually varying parameters, so a person cannot be held fixed ',
        'while its trajectory is redrawn: the factor applied to the initial ',
        'draw is triangular, and redrawing the process entries would move the ',
        "parameters too. Use method='moment' if the dynamics allow it.",
        call. = FALSE)
    }
    if (!.ctVarDecompCarrierIsDrawn(fit, layout, carrier, raw, base, augmented,
      nrows, design)) {
      stop("method='simulation' cannot give this model a between person ",
        "variance: this build's engine leaves every individually varying ",
        'parameter at its population value while generating states, so every ',
        'simulated person would be the same person and the between term would ',
        'come out at zero without anything having failed. Checked by ',
        'generating twice rather than assumed. Use method=\'moment\' if the ',
        'dynamics allow it.', call. = FALSE)
    }
  }

  nsubjects <- length(subjects)
  # With no carrier states every person is the same person, so redrawing them
  # buys nothing: one draw, and the paths are all this needs.
  draws <- if (length(carrier)) max(1L, ceiling(npersons / nsubjects)) else 1L

  nvar <- nmanifest + if (latents) nlatent else 0L
  npeople <- draws * nsubjects
  personmean <- matrix(NA_real_, npeople, nvar)
  persondet <- matrix(NA_real_, npeople, nvar)
  personstoch <- matrix(NA_real_, npeople, nvar)
  personmeas <- matrix(0, npeople, nvar)
  # The Monte Carlo corrections, accumulated per person and applied once.
  detcorrection <- matrix(0, npeople, nvar)
  betweencorrection <- matrix(0, npeople, nvar)

  at <- 0L
  for (drawi in seq_len(draws)) {
    pinned <- matrix(stats::rnorm(augmented * design$nsubjects), augmented,
      design$nsubjects)
    paths <- array(NA_real_, dim = c(npaths, nrows, augmented))
    for (path in seq_len(npaths)) {
      z <- stats::rnorm(layout$ndim)
      if (length(carrier)) for (si in seq_len(design$nsubjects)) {
        z[layout$zoffsets[si] + carrier] <- pinned[carrier, si]
      }
      drawn <- .ctBackendGenerateStates(fit, raw, z, base)
      paths[path, , ] <- t(matrix(as.numeric(drawn$states), augmented, nrows))
    }

    for (si in subjects) {
      at <- at + 1L
      rows <- which(design$subject == si)
      state <- as.numeric(population$T0MEANS)
      if (length(carrier)) state[carrier] <- paths[1L, rows[1L], carrier]
      mats <- .ctVarDecompPersonAt(fit,
        if (!is.null(tipreds)) as.numeric(tipreds[si, ]) else NULL, state, nlatent)
      # One npaths by length(rows) matrix per dynamic state. Built explicitly
      # rather than by indexing a three way array, so that a single path or a
      # single row cannot drop a dimension underneath the arithmetic.
      dynamic <- lapply(seq_len(nlatent), function(li)
        matrix(paths[, rows, li], npaths, length(rows)))

      for (vi in seq_len(nvar)) {
        latentonly <- vi > nmanifest
        g <- if (latentonly) dynamic[[vi - nmanifest]] else {
          linear <- matrix(mats$MANIFESTMEANS[vi, 1L], npaths, length(rows))
          for (li in seq_len(nlatent)) {
            linear <- linear + mats$LAMBDA[vi, li] * dynamic[[li]]
          }
          linear
        }
        conditional <- if (latentonly) 0 else
          .ctVarDecompSimulatedCondVar(g, type[vi], mats$MANIFESTcov[vi, vi], scale)
        if (!latentonly && !identical(scale, 'latent') && type[vi] != 0L) {
          g <- stats::plogis(g)
        }
        pathmean <- colMeans(g)
        overall <- rowMeans(g)
        centred <- g - overall
        personmean[at, vi] <- mean(overall)
        persondet[at, vi] <- .ctVarDecompPopVar(pathmean)
        personstoch[at, vi] <- mean(apply(g, 2L, stats::var))
        personmeas[at, vi] <- conditional
        detcorrection[at, vi] <- mean(apply(centred, 2L, stats::var)) / npaths
        betweencorrection[at, vi] <- stats::var(overall) / npaths
      }
    }
  }

  list(
    between = if (npeople > 1L)
      pmax(apply(personmean, 2L, stats::var) - colMeans(betweencorrection), 0)
      else rep(0, nvar),
    within.deterministic = pmax(colMeans(persondet) - colMeans(detcorrection), 0),
    within.stochastic = colMeans(personstoch),
    within.measurement = colMeans(personmeas),
    npersons = npeople)
}

# Does the engine's state generation actually draw the carrier states?
#
# It does, since `_ctsem_t0_factor!` -- before that the state pass built its
# initial factor from T0VAR alone, whose carrier entries are zero because the
# augmentation puts a random effect's spread in RAWPOPVAR, and so every
# generated subject got the population parameters. Generation is what this
# function asks about; the same factor is what a sampled `intoverstates=FALSE`
# fit draws its states through. This is the consumer side
# guard against that returning: the failure is a between person variance of
# zero with nothing having errored, which no amount of reading the
# specification would reveal. Two generation passes with the carrier entries of
# the initial draw far apart settle it.
#
# `test-julia-intoverstates.R` asserts the engine side invariant directly, that
# the state path's initial factor squares to the covariance the filter carries.
# This stays because it is what makes the refusal above honest at run time.
.ctVarDecompCarrierIsDrawn <- function(fit, layout, carrier, raw, base,
  augmented, nrows, design) {
  probe <- function(value) {
    z <- rep(0, layout$ndim)
    for (si in seq_len(design$nsubjects)) z[layout$zoffsets[si] + carrier] <- value
    states <- matrix(as.numeric(
      .ctBackendGenerateStates(fit, raw, z, base)$states), augmented, nrows)
    states[carrier, 1L]
  }
  any(abs(probe(5) - probe(-5)) > 1e-8)
}

# The measurement model's own variance at the drawn states. Nothing is
# integrated here: the linear predictor is known at each draw, so the
# conditional variance is read off it directly.
.ctVarDecompSimulatedCondVar <- function(linear, type, manifestvar, scale) {
  if (type == 0L) return(manifestvar)
  if (identical(scale, 'latent')) return(pi^2 / 3)
  p <- stats::plogis(linear)
  mean(p * (1 - p))
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
#' @param method How the model implied moments are obtained. \code{'moment'}
#'   runs the forward moment recursion, which is exact and needs the dynamics
#'   to be linear. \code{'simulation'} draws trajectories from the julia
#'   engine instead, which is what a model with state dependent \code{DRIFT}
#'   or \code{DIFFUSION} cells needs. \code{'auto'}, the default, is
#'   \code{'moment'} unless the model has such cells.
#' @param npaths Trajectories drawn per person when \code{method='simulation'}.
#'   Each person needs several: with one path, the average over that person's
#'   rows wanders with the path rather than estimating their mean, and the
#'   wander lands in the between person term. The two Monte Carlo corrections
#'   this makes necessary are applied, so the result is unbiased at any
#'   \code{npaths}, but small values are noisy.
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
#'   stationarity. The result is exact for a model with linear dynamics.
#'
#'   A model whose \code{DRIFT} or \code{DIFFUSION} cells depend on the state
#'   has no single transition to run that recursion with, so
#'   \code{method='simulation'} draws its trajectories from the engine and
#'   takes the same decomposition over the draws. It is refused rather than
#'   linearised at a point nothing chose. A state dependent \emph{measurement}
#'   cell -- \code{LAMBDA}, \code{MANIFESTMEANS}, \code{MANIFESTVAR} -- is
#'   refused by both routes, since turning a drawn state into an expected
#'   observation then needs the measurement model re-materialised at every row.
#'
#'   \strong{Random effect levels.} A model with a grouping level above the
#'   subject (\code{id = c('subject','study')} in \code{\link{ctModel}}) is
#'   fitted with \code{intoverpop='laplace'}, and gets one
#'   \code{between.<level>} column per level, each the variance of that
#'   level's own departure from the level outside it. They are orthogonal by
#'   construction and sum to \code{between}. Only the estimated (shrunk) route
#'   is available for such a fit, and the outer levels are shrunk hardest
#'   because they have the fewest units: on a 6-study, 30-subject example the
#'   study level came back at 0.24 against a generating 0.62 and the subject
#'   level at 0.03 against 0.22. Read those as lower bounds, and
#'   \code{summary(fit)} for the population parameters themselves.
#'
#'   \strong{Designs.} Nothing here assumes a balanced one. Each person is
#'   evaluated over their own observation times, so everyone measured at the
#'   same occasions, everyone at their own, and subjects with different numbers
#'   of occasions all work, as does an observation being missing -- a row still
#'   has a model implied distribution whether or not its indicator was seen, so
#'   every design row counts and the answer does not move with the missingness
#'   pattern. Persons are weighted equally rather than by row count: the
#'   estimand is over persons, so a frequently measured subject should not be a
#'   larger share of the population than a rarely measured one.
#'
#'   \code{method='simulation'} checks, by generating twice, that the engine
#'   draws a model's individually varying parameters rather than leaving every
#'   simulated person at the population values, and refuses rather than
#'   returning a between person variance of zero if it does not.
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
ctVarianceDecomposition <- function(fit, method = c('auto', 'moment', 'simulation'),
  persons = c('auto', 'model', 'estimated'), scale = c('latent', 'response'),
  npersons = 200L, npaths = 20L, latents = TRUE, quadpoints = 21L,
  subjects = 'all') {

  # `ctFit`, which both backends' fits carry, rather than naming the two
  # classes: this asks whether the argument is a fit at all, not which backend
  # produced it, and spelling it with a class literal puts it in front of the
  # duplication ratchet for a question it is not asking.
  if (!inherits(fit, 'ctFit')) {
    stop('fit object is not a ctsem fit!', call. = FALSE)
  }
  method <- match.arg(method)
  persons <- match.arg(persons)
  scale <- match.arg(scale)
  if (identical(method, 'auto')) {
    method <- if (nrow(.ctVarDecompStateDependent(fit))) 'simulation' else 'moment'
  }
  .ctVarDecompCheckModel(fit, method, scale)
  if (identical(method, 'simulation') && !.ctFitIsJulia(fit)) {
    stop("method='simulation' needs a backend='julia' fit: it draws the ",
      'trajectories from the engine, which the stan path has no entry point ',
      'for. Refit with backend=\'julia\'.', call. = FALSE)
  }
  if (identical(method, 'simulation') && !is.null(.ctBackendSpec(fit)$laplace)) {
    stop("method='simulation' is not available for a Laplace fit: its random ",
      'effects are not carrier states, so there is no part of the innovation ',
      'draw that pins a person while its path is redrawn.', call. = FALSE)
  }

  model <- .ctFitModelObject(fit)
  nlatent <- .ctFitNlatent(fit)
  continuoustime <- isTRUE(model$continuoustime)
  manifestNames <- model$manifestNames
  latentNames <- model$latentNames[seq_len(nlatent)]
  type <- as.integer(model$manifesttype)

  laplace <- .ctFitIsJulia(fit) && !is.null(.ctBackendSpec(fit)$laplace)
  if (identical(persons, 'auto')) {
    persons <- if (.ctFitIsJulia(fit) && !laplace) 'model' else 'estimated'
    if (laplace) {
      message("persons='estimated' for this Laplace fit: its random effects ",
        'are separate coordinates rather than carrier states, so a person ',
        'cannot be materialised at a drawn one. The between person variance ',
        'below is the spread of the estimated modes, which shrinkage ',
        'attenuates -- the outer levels hardest, since they have the fewest ',
        'units.')
    } else if (identical(persons, 'estimated')) {
      message("persons='estimated' for this stan fit: drawing persons from the ",
        'population distribution needs the model matrices materialised at a ',
        'drawn set of random effects, which only the julia engine can do. The ',
        'between person variance below is therefore the spread of the ',
        'estimated subjects, which shrinkage attenuates.')
    }
  }
  if (identical(persons, 'model') && laplace) {
    stop("persons='model' is not available for a Laplace fit: drawing a person ",
      "needs a draw from each level's own random effect covariance, and the ",
      'engine materialises matrices at a carrier state, which a Laplace fit ',
      "has none of. Use persons='estimated', reading its between person ",
      'variance as attenuated by shrinkage.', call. = FALSE)
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

  nmanifest <- length(manifestNames)
  components <- if (identical(method, 'simulation')) {
    if (identical(persons, 'estimated')) {
      stop("method='simulation' draws its persons from the population ",
        "distribution, so persons='estimated' has nothing to mean here: a ",
        "subject's estimated random effects would have to be turned back into ",
        'the standard normals the engine draws them from, through a factor ',
        "this does not have. Use persons='model'.", call. = FALSE)
    }
    .ctVarDecompSimulation(fit, design, nlatent, latents, type, scale,
      as.integer(npersons), as.integer(npaths), nmanifest, wanted)
  } else if (laplace) {
    gh <- if (identical(scale, 'response')) .ctVarDecompGaussHermite(quadpoints) else NULL
    .ctVarDecompLaplaceComponents(fit, design, nlatent, latents, type, scale,
      gh, continuoustime, nmanifest, wanted)
  } else {
    people <- if (.ctFitIsJulia(fit)) {
      .ctVarDecompJuliaPersons(fit, design, nlatent, persons,
        as.integer(npersons), wanted)
    } else .ctVarDecompStanPersons(fit, design, nlatent, wanted)
    people <- people[vapply(people, function(p) length(p$rows) > 0L, logical(1L))]
    if (!length(people)) stop('No usable persons.', call. = FALSE)
    gh <- if (identical(scale, 'response')) .ctVarDecompGaussHermite(quadpoints) else NULL
    .ctVarDecompMomentComponents(fit, design, people, nlatent, latents, type,
      scale, gh, continuoustime, nmanifest)
  }

  out <- data.frame(
    variable = c(manifestNames, if (latents) latentNames),
    type = c(rep('manifest', nmanifest), if (latents) rep('latent', nlatent)),
    between = components$between,
    within.deterministic = components$within.deterministic,
    within.stochastic = components$within.stochastic,
    within.measurement = components$within.measurement,
    stringsAsFactors = FALSE)
  # One column per random effect level, where there is more than one: with a
  # single level it would repeat `between` under another name.
  if (!is.null(components$levels) && nrow(components$levels) > 1L) {
    for (l in rownames(components$levels)) {
      out[[paste0('between.', l)]] <- components$levels[l, ]
    }
  }
  out$within <- out$within.deterministic + out$within.stochastic + out$within.measurement
  out$total <- out$between + out$within
  # A variable with no variance at all has no proportions, and NA says that
  # where 0/0 would print NaN. A constant auxiliary state is the case: a
  # higher order model's carried coordinate contributes nothing to anything.
  for (part in c('between', 'within.deterministic', 'within.stochastic',
    'within.measurement', 'within')) {
    out[[paste0('prop.', part)]] <- ifelse(out$total > 0,
      out[[part]] / out$total, NA_real_)
  }
  rownames(out) <- NULL

  attr(out, 'method') <- method
  attr(out, 'persons') <- persons
  attr(out, 'scale') <- scale
  attr(out, 'npersons') <- components$npersons
  attr(out, 'npaths') <- if (identical(method, 'simulation')) as.integer(npaths)
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
  levelcols <- grep('^between\\.', names(x), value = TRUE)
  if (length(levelcols)) {
    cat('\nBetween person variance by level\n')
    print(data.frame(variable = x$variable,
      as.data.frame(lapply(as.data.frame(x)[levelcols], round, digits))),
      row.names = FALSE)
  }
  cat('\nComputed by the ', attr(x, 'method'), ' route',
    if (identical(attr(x, 'method'), 'simulation'))
      paste0(', ', attr(x, 'npaths'), ' paths per person') else '', '.\n', sep = '')
  cat('Between person variance refers to ',
    if (identical(attr(x, 'persons'), 'model'))
      paste0(attr(x, 'npersons'), ' persons drawn from the fitted population distribution')
    else paste0('the ', attr(x, 'npersons'),
      ' estimated subjects, attenuated by shrinkage'),
    '.\n', sep = '')
  cat('Within person variance is over each person\'s own observation times; ',
    'the deterministic part is the mean path moving.\n', sep = '')
  if (identical(attr(x, 'method'), 'simulation')) {
    cat('The dynamics are integrated along drawn trajectories rather than ',
      'evaluated at one state, so no evaluation point is reported; the ',
      'measurement model does not depend on the state, which was checked.\n',
      sep = '')
  }
  if (identical(attr(x, 'scale'), 'latent')) {
    cat("Non-Gaussian indicators are decomposed on the latent response scale ",
      '(measurement variance pi^2/3).\n', sep = '')
  } else cat('Non-Gaussian indicators are decomposed on the response scale.\n')
  invisible(x)
}
