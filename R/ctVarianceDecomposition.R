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
# Uncertainty. The decomposition is a function of the parameters, so it is
# computed at each of a set of the fit's parameter draws (normal approximation
# draws on an optimised fit, posterior draws on a sampled one), and each
# quantity is reported as the median of its own draws, with its 2.5% and 97.5%
# quantiles beside it. A quantity built from other quantities' medians would
# not be the median of anything, so the parts' medians need not add up to the
# total's. Every person at every draw runs the same recursion on matrices of
# the same size, so it runs over all of them at once (the batched arithmetic
# below), and the persons are the same at every draw: the error of drawing them
# is then nearly the same at each draw, and moves the draws together rather
# than widening their spread. Before the batching, one call at the estimate
# took 10 s on the book's RI-CLPM (dev3), nearly all of it one engine call per
# person for their initial means.
#
# Both backends are served from one body. Everything read here -- the
# parameter matrices from the engine, stan's saved subject matrices, the design
# accessors in R/ctBackendKalman.R, `.ctFitConditionalCells()` -- already means
# the same thing on each.


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


# Batched small-matrix arithmetic ---------------------------------------------
#
# An array x[k, i, j] holds K matrices -- one per (draw, person) trajectory, or
# one per distinct interval -- and each operation below is a few vector
# operations of length K times the matrix size. Every trajectory runs the same
# recursion on matrices of the same size, so running them together costs R's
# per-call overhead once rather than once a person and draw. Nothing in base R
# or the packages ctsem depends on multiplies stacks of matrices.

# a: K x p x q, b: K x q x r; K x p x r.
.ctVarDecompMul <- function(a, b) {
  da <- dim(a)
  K <- da[1L]; p <- da[2L]; q <- da[3L]; r <- dim(b)[3L]
  # Entry (k, i, l) of the result reads b[k, j, l] for every i.
  index <- rep(seq_len(K), p * r) + K * rep(seq_len(r) - 1L, each = K * p)
  out <- numeric(K * p * r)
  for (j in seq_len(q)) {
    out <- out + as.vector(a[, , j]) * as.vector(b[, j, ])[index]
  }
  array(out, c(K, p, r))
}

# a: K x p x q, x: K x q; K x p.
.ctVarDecompMulVec <- function(a, x) {
  d <- dim(a)
  out <- numeric(d[1L] * d[2L])
  for (j in seq_len(d[3L])) out <- out + as.vector(a[, , j]) * x[, j]
  matrix(out, d[1L], d[2L])
}

.ctVarDecompTr <- function(a) aperm(a, c(1L, 3L, 2L))

.ctVarDecompSym <- function(a) (a + .ctVarDecompTr(a)) / 2

# A P A', symmetrised as the unbatched recursion did after every row.
.ctVarDecompSandwich <- function(a, p) {
  .ctVarDecompSym(.ctVarDecompMul(.ctVarDecompMul(a, p), .ctVarDecompTr(a)))
}

.ctVarDecompEye <- function(K, n) array(rep(as.vector(diag(n)), each = K), c(K, n, n))

# Each matrix's diagonal, K x n.
.ctVarDecompDiag <- function(a) {
  d <- dim(a)
  i <- rep(seq_len(d[2L]), each = d[1L])
  matrix(a[cbind(rep(seq_len(d[1L]), d[2L]), i, i)], d[1L], d[2L])
}

# diag(L P L'), K x m.
.ctVarDecompQuadDiag <- function(l, p) {
  matrix(rowSums(.ctVarDecompMul(l, p) * l, dims = 2L), dim(l)[1L], dim(l)[2L])
}

# Each matrix's 1-norm.
.ctVarDecompNorm <- function(a) {
  sums <- matrix(rowSums(aperm(abs(a), c(1L, 3L, 2L)), dims = 2L), dim(a)[1L])
  do.call(pmax, as.data.frame(sums))
}

# The exponential of matrices whose 1-norm is at most a half, by Taylor
# series: the remainder after degree 14 is under 0.5^15 / 15!, about 2e-17.
# `.ctVarDecompSteps()` scales every interval down to there first.
.ctVarDecompExpmSmall <- function(x, degree = 14L) {
  eye <- .ctVarDecompEye(dim(x)[1L], dim(x)[2L])
  out <- eye
  for (j in rev(seq_len(degree))) out <- eye + .ctVarDecompMul(x, out) / j
  out
}


# Each interval's transition, intercept and innovation covariances ------------
#
# For K intervals at once: `drift` K x n x n, `cint` K x n, `diffusions` a list
# of K x n x n (the whole DIFFUSIONcov, and with `sources` each source's share
# of it), `dt` length K. Returns the transition K x n x n, the intercept K x n
# and one innovation covariance per diffusion.
#
# Continuous time is scaling and squaring of the moments themselves. Over an
# interval short enough that every block below has norm at most a half, the
# block exponentials give the step to rounding: expm([[A, c], [0, 0]] t), whose
# top right block is int_0^t expm(A s) ds c and so is defined when DRIFT is
# singular, where asymCINT is not; and Van Loan's expm([[A, Q], [0, -A']] t)
# for the innovation, as `.ctNetworkInnovation()` (R/ctGraph.R) computes it,
# needing no stationarity assumption. An interval twice as long is two of them
# in a row,
#
#   A <- A A,   b <- A b + b,   Q <- A Q A' + Q,
#
# which is the moment recursion over two equal intervals, so the doubling
# carries only the process's own decay. Exponentiating the Van Loan block over
# the whole interval instead carries expm(-A' dt), which grows as fast as the
# process decays.
#
# Discrete time composes the one-step moments the interval's number of times.
.ctVarDecompSteps <- function(drift, cint, diffusions, dt, continuoustime) {
  K <- dim(drift)[1L]
  n <- dim(drift)[2L]
  latent <- seq_len(n)
  if (!continuoustime) {
    steps <- as.integer(round(dt))
    if (any(steps < 1L)) stop(call. = FALSE,
      'A discrete time model needs whole number time steps; found an interval of ',
      dt[steps < 1L][1L], '.')
    transition <- .ctVarDecompEye(K, n)
    intercept <- matrix(0, K, n)
    innovation <- lapply(diffusions, function(q) array(0, dim(q)))
    for (j in seq_len(max(steps))) {
      i <- which(steps >= j)
      a <- drift[i, , , drop = FALSE]
      for (g in seq_along(innovation)) {
        innovation[[g]][i, , ] <- .ctVarDecompSandwich(a,
          innovation[[g]][i, , , drop = FALSE]) + diffusions[[g]][i, , , drop = FALSE]
      }
      intercept[i, ] <- .ctVarDecompMulVec(a, intercept[i, , drop = FALSE]) +
        cint[i, , drop = FALSE]
      transition[i, , ] <- .ctVarDecompMul(a, transition[i, , , drop = FALSE])
    }
    return(list(transition = transition, intercept = intercept,
      innovation = innovation))
  }

  # A bound on every block's norm, and the halvings that bring it under a half.
  size <- dt * (.ctVarDecompNorm(drift) + .ctVarDecompNorm(.ctVarDecompTr(drift)) +
    rowSums(abs(cint)) + do.call(pmax, c(lapply(diffusions, .ctVarDecompNorm), 0)))
  halvings <- pmax(0, ceiling(log2(size / 0.5)))
  halvings[!is.finite(halvings)] <- 0
  short <- dt / 2^halvings

  block <- array(0, c(K, n + 1L, n + 1L))
  block[, latent, latent] <- drift * short
  block[, latent, n + 1L] <- cint * short
  e <- .ctVarDecompExpmSmall(block)
  transition <- e[, latent, latent, drop = FALSE]
  intercept <- matrix(e[, latent, n + 1L], K, n)
  innovation <- lapply(diffusions, function(q) {
    vanloan <- array(0, c(K, 2L * n, 2L * n))
    vanloan[, latent, latent] <- drift * short
    vanloan[, latent, n + latent] <- q * short
    vanloan[, n + latent, n + latent] <- -.ctVarDecompTr(drift) * short
    e <- .ctVarDecompExpmSmall(vanloan)
    .ctVarDecompSym(.ctVarDecompMul(e[, latent, n + latent, drop = FALSE],
      .ctVarDecompTr(e[, latent, latent, drop = FALSE])))
  })
  for (j in seq_len(max(halvings, 0))) {
    i <- which(halvings >= j)
    a <- transition[i, , , drop = FALSE]
    for (g in seq_along(innovation)) {
      before <- innovation[[g]][i, , , drop = FALSE]
      innovation[[g]][i, , ] <- .ctVarDecompSandwich(a, before) + before
    }
    intercept[i, ] <- .ctVarDecompMulVec(a, intercept[i, , drop = FALSE]) +
      intercept[i, , drop = FALSE]
    transition[i, , ] <- .ctVarDecompMul(a, a)
  }
  list(transition = transition, intercept = intercept, innovation = innovation)
}


# Which trajectories share their dynamics -------------------------------------
#
# One group per distinct row of `key` (each trajectory's DRIFT, DIFFUSIONcov
# and CINT), so an interval's step is computed once for every trajectory that
# has the same matrices: once a draw where only intercepts vary between
# persons, once a person and draw where the dynamics do. Rows are matched by a
# projection, checked exactly, and by their text where the check fails.
.ctVarDecompGroups <- function(key) {
  if (!nrow(key)) return(integer())
  hash <- as.vector(key %*% (1 / sqrt(seq_len(ncol(key)) + 0.37)))
  group <- match(hash, unique(hash))
  first <- match(seq_len(max(group)), group)
  if (isTRUE(all(key == key[first[group], , drop = FALSE]))) return(group)
  text <- do.call(paste, c(as.data.frame(key), sep = '\r'))
  match(text, unique(text))
}


# The marginal moments of every trajectory over its rows -----------------------
#
# No data enters: these are the model's moments given that person's parameters,
# which is what a decomposition of the model's implied variance needs. A
# smoothed trajectory would answer a different question and would be shrunk
# toward the observations.
#
# A trajectory is one person at one parameter draw. `mats` holds each one's
# matrices, K x r x c arrays with the draw running fastest (trajectory
# d + ndraws (p - 1)), and `times` and `tdpreds` each person's rows. Row r of
# every trajectory with at least r rows advances together, and each step is
# computed once for each distinct pair of dynamics group and interval.
#
# `covariance = FALSE` carries the means alone, which is all the outer levels
# need on the latent scale. With `sources`, the state covariance is also kept
# apart by where it came from: the initial state, and each group of latents'
# innovations (`.ctVarDecompSources`). They add up because the groups'
# innovations are independent and the transition is linear. Beside each, the
# sum over pairs of rows of the state's covariance, `K_t = P_t + A_t K_{t-1}`
# being a row's covariance with every row before it and itself, which is what
# the variance of the record's mean is made of.
.ctVarDecompMoments <- function(mats, ndraws, times, tdpreds, continuoustime,
  covariance = TRUE, sources = NULL) {
  B <- dim(mats$DRIFT)[1L]
  n <- dim(mats$DRIFT)[2L]
  rows <- lengths(times)
  owner <- rep(seq_along(times), each = ndraws)
  count <- rows[owner]
  maxrows <- max(rows)
  active <- lapply(seq_len(maxrows), function(r) which(count >= r))
  grid <- matrix(NA_real_, length(times), maxrows)
  grid[cbind(rep(seq_along(times), rows), sequence(rows))] <- unlist(times)
  split <- length(sources) > 0L

  # Every interval each trajectory steps over, by dynamics group.
  group <- .ctVarDecompGroups(cbind(matrix(mats$DRIFT, B),
    matrix(mats$DIFFUSIONcov, B), matrix(mats$CINT, B)))
  ngroup <- max(group, 1L)
  later <- seq_len(maxrows)[-1L]
  gaps <- lapply(later, function(r) {
    b <- active[[r]]
    grid[owner[b], r] - grid[owner[b], r - 1L]
  })
  intervals <- unique(unlist(gaps))
  keys <- lapply(seq_along(later), function(i)
    group[active[[later[i]]]] + ngroup * (match(gaps[[i]], intervals) - 1))
  if (length(later)) {
    allkeys <- unlist(keys)
    table <- unique(allkeys)
    first <- unlist(active[later])[match(table, allkeys)]
    diffusions <- if (covariance || split) {
      c(list(mats$DIFFUSIONcov[first, , , drop = FALSE]),
        if (split) lapply(sources, function(s) {
          d <- array(0, c(length(first), n, n))
          d[, s, s] <- mats$DIFFUSIONcov[first, s, s, drop = FALSE]
          d
        }))
    } else list()
    steps <- .ctVarDecompSteps(mats$DRIFT[first, , , drop = FALSE],
      matrix(mats$CINT[first, , , drop = FALSE], length(first), n), diffusions,
      intervals[(table - 1) %/% ngroup + 1], continuoustime)
  }

  state <- matrix(mats$T0MEANS, B, n)
  cov <- if (covariance) .ctVarDecompSym(mats$T0cov)
  if (split) {
    part <- c(list(.ctVarDecompSym(mats$T0cov)),
      lapply(unname(sources), function(s) array(0, c(B, n, n))))
    carried <- lapply(part, function(x) array(0, c(B, n, n)))
    pairs <- carried
  }
  shifted <- !is.null(tdpreds) && !is.null(mats$TDPREDEFFECT)
  out <- vector('list', maxrows)
  for (r in seq_len(maxrows)) {
    b <- active[[r]]
    if (r > 1L) {
      k <- match(keys[[r - 1L]], table)
      a <- steps$transition[k, , , drop = FALSE]
      state[b, ] <- .ctVarDecompMulVec(a, state[b, , drop = FALSE]) +
        steps$intercept[k, , drop = FALSE]
      if (covariance) cov[b, , ] <- .ctVarDecompSandwich(a,
        cov[b, , , drop = FALSE]) + steps$innovation[[1L]][k, , , drop = FALSE]
      if (split) for (g in seq_along(part)) {
        part[[g]][b, , ] <- .ctVarDecompSandwich(a, part[[g]][b, , , drop = FALSE]) +
          if (g > 1L) steps$innovation[[g]][k, , , drop = FALSE] else 0
        carried[[g]][b, , ] <- .ctVarDecompMul(a, carried[[g]][b, , , drop = FALSE])
      }
    }
    if (split) for (g in seq_along(part)) {
      carried[[g]][b, , ] <- carried[[g]][b, , , drop = FALSE] +
        part[[g]][b, , , drop = FALSE]
      pairs[[g]][b, , ] <- pairs[[g]][b, , , drop = FALSE] +
        carried[[g]][b, , , drop = FALSE] +
        .ctVarDecompTr(carried[[g]][b, , , drop = FALSE]) - part[[g]][b, , , drop = FALSE]
    }
    # A time dependent predictor is an instantaneous shift at the observed row,
    # applied after the time update -- matching the filter, which does
    # `state[1:nlatent] += (TDPREDEFFECT * tdpreds[rowi])'` at that point.
    if (shifted) {
      persons <- unique(owner[b])
      values <- matrix(vapply(persons, function(p) tdpreds[[p]][r, ],
        numeric(ncol(tdpreds[[1L]]))), ncol = ncol(tdpreds[[1L]]), byrow = TRUE)
      state[b, ] <- state[b, , drop = FALSE] + .ctVarDecompMulVec(
        mats$TDPREDEFFECT[b, , , drop = FALSE],
        values[rep(seq_along(persons), each = ndraws), , drop = FALSE])
    }
    lambda <- mats$LAMBDA[b, , , drop = FALSE]
    out[[r]] <- list(index = b, latentmean = state[b, , drop = FALSE],
      linearmean = .ctVarDecompMulVec(lambda, state[b, , drop = FALSE]) +
        matrix(mats$MANIFESTMEANS[b, , , drop = FALSE], length(b)))
    if (covariance) {
      out[[r]]$latentvar <- .ctVarDecompDiag(cov[b, , , drop = FALSE])
      out[[r]]$linearvar <- .ctVarDecompQuadDiag(lambda, cov[b, , , drop = FALSE])
    }
    if (split) {
      out[[r]]$partlatent <- vapply(part, function(x)
        .ctVarDecompDiag(x[b, , , drop = FALSE]), matrix(0, length(b), n))
      out[[r]]$partlinear <- vapply(part, function(x) .ctVarDecompQuadDiag(lambda,
        x[b, , , drop = FALSE]), matrix(0, length(b), dim(lambda)[2L]))
    }
  }
  result <- list(rows = out)
  if (split) {
    # The record mean's variance from each part: the sum over pairs of rows
    # over the number of pairs.
    result$levellatent <- vapply(pairs, function(p) .ctVarDecompDiag(p) / count^2,
      matrix(0, B, n))
    result$levellinear <- vapply(pairs, function(p)
      .ctVarDecompQuadDiag(mats$LAMBDA, p) / count^2,
      matrix(0, B, dim(mats$LAMBDA)[2L]))
  }
  result
}

# `sources` resolved against the model: a named list of latent names, each a
# group whose innovations are independent of every other group's, checked on
# the fitted DIFFUSIONcov. Returned as latent indices.
.ctVarDecompSources <- function(fit, sources, nlatent) {
  if (is.null(sources)) return(NULL)
  latents <- .ctFitModelObject(fit)$latentNames[seq_len(nlatent)]
  if (!is.list(sources) || is.null(names(sources)) || any(!nzchar(names(sources))))
    stop("sources must be a named list of latent names, one entry per group of ",
      "processes whose innovations are independent of the others'.", call. = FALSE)
  if (any(names(sources) %in% c("initial", "level"))) stop("'initial' and ",
    "'level' name parts the decomposition reports; call the sources otherwise.",
    call. = FALSE)
  index <- lapply(sources, function(s) {
    i <- match(s, latents)
    if (anyNA(i)) stop("sources names latents this model does not have: ",
      paste(s[is.na(i)], collapse = ", "), ".", call. = FALSE)
    i
  })
  all <- unlist(index)
  if (anyDuplicated(all)) stop("A latent appears in more than one source.",
    call. = FALSE)
  diffusion <- as.matrix(suppressMessages(ctSummaryMatrices(fit))$DIFFUSIONcov)
  outside <- setdiff(seq_len(nlatent), all)
  if (length(outside) && any(abs(diffusion[outside, ]) > 0)) stop(
    "Latents with innovations of their own are in no source: ",
    paste(latents[outside], collapse = ", "), ".", call. = FALSE)
  for (a in seq_along(index)) for (b in seq_along(index)) {
    if (a == b) next
    if (any(abs(diffusion[index[[a]], index[[b]]]) >
        1e-8 * max(abs(diffusion)))) stop("The innovations of sources '",
      names(sources)[a], "' and '", names(sources)[b], "' are correlated in ",
      "the fitted DIFFUSIONcov, so the variance does not split between them. ",
      "Put the correlated latents in one source.", call. = FALSE)
  }
  index
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
# measurement variance: MANIFESTcov for a Gaussian indicator (one value per
# trajectory, or one for all), and pi^2/3 for a binary or ordinal one, the
# logistic variance implied by the threshold form
# P(y <= k | eta) = inv_logit(tau_k - eta).
#
# On the response scale the link is integrated over the row's own normal
# marginal by Gauss-Hermite, which is exact to quadrature error.
.ctVarDecompMeasurement <- function(mean, var, type, manifestvar, scale, gh) {
  if (type == 0L) {
    return(list(expected = mean, varmean = var,
      condvar = rep_len(manifestvar, length(mean))))
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


# Parameter draws -------------------------------------------------------------
#
# The draws each quantity is computed at: `nsamples` of the fit's own, evenly
# spaced as ctExtract(nsamples = ) takes them, or the estimate alone when
# intervals are turned off, when the fit has no draws, or when `pointonly`
# gives the reason a route cannot use them.
.ctVarDecompDraws <- function(fit, nsamples, pointonly = NULL) {
  if (length(nsamples) != 1L || is.na(nsamples) || nsamples < 0) {
    stop('nsamples must be 0, for the estimate without intervals, or a number ',
      'of parameter draws.', call. = FALSE)
  }
  julia <- .ctFitIsJulia(fit)
  point <- list(raw = matrix(as.numeric(if (julia) fit$estimate$raw else
    fit$stanfit$rawest), 1L), rows = NULL, intervals = FALSE)
  if (nsamples < 1) return(point)
  if (!is.null(pointonly)) {
    message(pointonly)
    return(point)
  }
  samples <- if (julia) fit$estimate$rawposterior else fit$stanfit$rawposterior
  if (is.null(samples) || length(dim(samples)) != 2L || nrow(samples) < 2L) {
    message('This fit has no parameter draws, so values are at the estimate, ',
      'without intervals.')
    return(point)
  }
  rows <- unique(round(seq(1, nrow(samples), length.out =
    min(nrow(samples), as.integer(nsamples)))))
  list(raw = matrix(as.numeric(samples[rows, , drop = FALSE]), length(rows)),
    rows = rows, intervals = TRUE)
}


# Persons ---------------------------------------------------------------------
#
# Who the persons are is decided once, with every random draw it needs
# (`.ctVarDecompPlan()`), and their matrices are built at each chunk of
# parameter draws (`.ctVarDecompMaterialise()`). A chunk is a batch:
#
#   rows    each person's design rows
#   unit    persons x levels: the unit of each random effect level a person
#           belongs to
#   mats    one entry per step, each a named list of arrays over the
#           trajectories, draw fastest: the matrices with the effects of that
#           level and every level outside it, innermost first, the last entry
#           being the population's
#   ndraws  the draws in the chunk
#
# A level's own contribution is then the difference between successive entries
# of `mats`, which is what makes the between person variance a sum over levels
# and what makes the three representations interchangeable here. It also takes
# the design out of the between term: two persons measured at different
# occasions have different time-averaged expected values even with identical
# parameters, and subtracting the population value at that same design removes
# it rather than reporting it as an individual difference.
#
# Each subject's design is used `round(npersons / nsubjects)` times, at least
# once, so every design counts equally.

# Standard normal draws, n x k, standardised to exactly zero mean and identity
# covariance over the rows where there are enough rows: the persons drawn then
# have the population's mean and covariance exactly, so a between person
# variance that is linear in the effects with the same weights for everyone
# (random intercepts) has no sampling error, and others less.
.ctVarDecompNormals <- function(n, k) {
  z <- matrix(stats::rnorm(n * k), n, k)
  if (k < 1L || n < k + 2L) return(z)
  z <- sweep(z, 2L, colMeans(z))
  z %*% solve(chol(crossprod(z) / (n - 1L)))
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

# Trim the augmented state out of matrices, each an array with one leading
# entry per trajectory.
#
# The subject and population matrices are already over the real latent
# processes -- except the T0 pair, which keeps the carrier rows because that is
# where the random effects live -- so those are cut here. Matrices
# materialised at a carrier state are untrimmed, so that the carrier block can
# be read whole, and are cut by the same rule .ctBackendTrimAugmented()
# applies.
.ctVarDecompTrim <- function(mats, nlatent, augmented) {
  keep <- seq_len(nlatent)
  if (augmented) {
    for (name in intersect(names(mats), c('DRIFT', 'DIFFUSIONcov'))) {
      mats[[name]] <- mats[[name]][, keep, keep, drop = FALSE]
    }
    for (name in intersect(names(mats), c('CINT', 'TDPREDEFFECT'))) {
      mats[[name]] <- mats[[name]][, keep, , drop = FALSE]
    }
    if (!is.null(mats$LAMBDA)) mats$LAMBDA <- mats$LAMBDA[, , keep, drop = FALSE]
  }
  mats$T0MEANS <- mats$T0MEANS[, keep, , drop = FALSE]
  mats$T0cov <- mats$T0cov[, keep, keep, drop = FALSE]
  mats
}

# The flat positions of the named matrices in the engine's layout, and the
# engine's columns at those positions as one array per matrix: columns x rows
# x columns of the matrix.
.ctVarDecompFlatRows <- function(layout, names) {
  unlist(lapply(names, function(name) {
    i <- match(name, layout$matrix)
    layout$offset[i] + seq_len(layout$nrow[i] * layout$ncol[i])
  }))
}

.ctVarDecompFlatArrays <- function(flat, layout, names) {
  out <- list()
  at <- 0L
  for (name in names) {
    i <- match(name, layout$matrix)
    size <- layout$nrow[i] * layout$ncol[i]
    out[[name]] <- array(t(flat[at + seq_len(size), , drop = FALSE]),
      c(ncol(flat), layout$nrow[i], layout$ncol[i]))
    at <- at + size
  }
  out
}

# Arrays for some trajectories written into arrays over all of them.
.ctVarDecompScatter <- function(into, at, mats, total) {
  for (name in names(mats)) {
    if (is.null(into[[name]])) into[[name]] <- array(NA_real_,
      c(total, dim(mats[[name]])[-1L]))
    into[[name]][at, , ] <- mats[[name]]
  }
  into
}

# The model matrices at raw vectors and whole filter states, a column each.
# `.ctBackendParMatricesFlat()` takes either a state of the latent processes a
# column, or one whole filter state for every column; a carrier person needs
# their own carrier values at every draw.
.ctVarDecompAtStates <- function(fit, raw, states, tipreds, rows) {
  storage.mode(raw) <- 'double'
  storage.mode(states) <- 'double'
  module <- .ctJuliaModule(.ctBackendSpec(fit)$project)
  arguments <- list(.ctJuliaObjective(fit), .ctJuliaPut(raw), time = 0, dt = 0,
    state = .ctJuliaPut(states), rows = .ctJuliaVector(as.integer(rows)))
  if (length(tipreds)) arguments$tipreds <- .ctJuliaVector(as.numeric(tipreds))
  result <- .ctBackendJuliaValue(do.call(module$ctsem_parameter_matrices, arguments))
  matrix(as.numeric(result), nrow = length(rows))
}

# Who the persons are, and every random draw that makes them. The persons are
# each subject's design `replicates` times over, and they and their draws are
# the same at every parameter draw.
.ctVarDecompPlan <- function(fit, design, nlatent, persons, npersons, subjects,
  levels, drawrows = NULL) {
  julia <- .ctFitIsJulia(fit)
  nlevels <- length(levels)
  coordinates <- julia && .ctSpecEffectsAreCoordinates(.ctBackendSpec(fit))
  route <- if (!julia) 'stan' else if (coordinates) 'coordinate' else 'carrier'
  plan <- list(route = route, persons = persons, nlatent = nlatent,
    nlevels = nlevels)

  if (identical(route, 'coordinate')) {
    # A joint sample that kept its draws of the effects pairs each parameter
    # draw with its own, as `.ctBackendPlausibleValues()` does. At the
    # estimate there is no draw to pair.
    saved <- fit$sample$effects
    plan$saved <- identical(persons, 'estimated') && length(drawrows) > 0L &&
      identical(fit$sample$target, 'joint') && length(fit$sample$effect_mean) > 0L &&
      !is.null(saved) && nrow(saved) == nrow(.ctBackendRawSamples(fit))
  }
  replicates <- if (!julia || isTRUE(plan$saved)) 1L else
    max(1L, as.integer(round(npersons / length(subjects))))
  plan$replicates <- replicates
  plan$subjects <- rep(subjects, replicates)
  copy <- rep(seq_len(replicates), each = length(subjects))
  plan$copy <- copy
  plan$rows <- lapply(plan$subjects, function(si) which(design$subject == si))
  npeople <- length(plan$subjects)

  # Time independent predictors, a call to the engine for each distinct row.
  tipreds <- if (julia && length(.ctFitModelObject(fit)$TIpredNames))
    .ctFitTIpredData(fit) else NULL
  if (is.null(tipreds)) {
    plan$group <- rep(1L, npeople)
    plan$tipreds <- list(NULL)
  } else {
    key <- do.call(paste, c(as.data.frame(tipreds[plan$subjects, , drop = FALSE]),
      sep = '\r'))
    plan$group <- match(key, unique(key))
    plan$tipreds <- lapply(match(seq_len(max(plan$group)), plan$group),
      function(p) as.numeric(tipreds[plan$subjects[p], ]))
  }

  if (identical(route, 'coordinate')) {
    structure <- .ctFitRandomEffectLevels(fit)
    plan$unit <- matrix(vapply(seq_len(nlevels), function(l)
      as.integer(structure[[l]]$units[plan$subjects] +
        (copy - 1L) * structure[[l]]$nunits), integer(npeople)), npeople, nlevels)
    plan$position <- .ctVarDecompPositions(fit, levels)
    if (identical(persons, 'model')) plan$normals <- lapply(seq_len(nlevels),
      function(l) .ctVarDecompNormals(structure[[l]]$nunits * replicates,
        length(plan$position[[l]])))
    return(plan)
  }
  # One level at most here, whose units are the persons.
  plan$unit <- matrix(rep(seq_len(npeople), nlevels), npeople, nlevels)
  if (identical(route, 'stan')) return(plan)

  # An individually varying parameter is carried as a latent state with no
  # drift and no diffusion, so a person *is* a carrier vector and the engine
  # materialises their matrices at it, applying every transform through the
  # same code the likelihood uses. persons='model' draws the vector from the
  # population distribution (the carrier block of the augmented T0MEANS and
  # T0cov, shifted by that person's time independent predictors);
  # persons='estimated' reads each subject's own carrier values off the filter
  # instead, at plausible values: the smoothed value plus a draw from the
  # posterior the fit stored (`.ctVarDecompCarrierNoise()`), so the spread
  # over persons is not shrunk.
  carrier <- .ctVarDecompCarrier(fit, nlatent)
  # The model says it has individual differences and this route cannot see
  # them. Refused rather than returned, because what it would return is a
  # between person variance of zero, which is what a model with no individual
  # differences correctly returns -- so nothing downstream could tell the two
  # apart.
  if (!length(carrier) && nlevels) {
    stop('This fit declares individually varying parameters and they are ',
      'neither carrier states nor separate coordinates, so no route here can ',
      'see them and the between person variance would come out at zero. This ',
      'is a bug rather than a limitation of the model -- please report the fit.',
      call. = FALSE)
  }
  plan$carrier <- carrier
  if (length(carrier)) {
    if (identical(persons, 'model')) {
      plan$normals <- .ctVarDecompNormals(npeople, length(carrier))
    } else {
      smoothed <- .ctVarDecompSubjectMatrices(fit)$subj_T0MEANS
      plan$states <- matrix(smoothed[1L, plan$subjects, carrier, 1L], npeople) +
        .ctVarDecompCarrierNoise(fit, plan$subjects, carrier)
    }
  }
  plan
}

# Each person's matrices at a chunk of parameter draws (`raw`, draws x
# parameters, the fit's rows `drawrows` of its draws, or the estimate).
.ctVarDecompMaterialise <- function(plan, fit, raw, drawrows = NULL) {
  switch(plan$route,
    stan = .ctVarDecompStanMats(plan, fit, raw),
    carrier = .ctVarDecompCarrierMats(plan, fit, raw),
    coordinate = .ctVarDecompCoordinateMats(plan, fit, raw, drawrows))
}

# The batch from per-step matrices.
.ctVarDecompBatch <- function(plan, mats, ndraws) {
  list(rows = plan$rows, unit = plan$unit, ndraws = ndraws, mats = mats)
}

# The needed matrices this model's engine layout has.
.ctVarDecompLayoutNeeds <- function(layout) {
  intersect(.ctVarDecompNeeded, layout$matrix)
}


# Coordinate fits: a person is an effect vector --------------------------------
#
# `intoverpop='laplace'`, `'variational'` and `'none'` keep the random effects
# as coordinates rather than states, so there is no carrier to materialise at.
# A person is instead a raw parameter vector per level of effects: the
# population vector shifted by the effects of that level and every level
# outside it, at the effects' positions in the raw vector, which the engine
# then materialises with the subject's time independent predictors. The vector
# is the one before TI-predictor effects, which the engine adds from the
# subject's own predictors.
#
# persons='model' draws each unit's effects from its level's population
# covariance at each parameter draw, from standard normals drawn once:
# `root %*% z`, the root symmetric (`.ctSymRoot()`) so that it moves smoothly
# with the parameters. The draw is per *unit*: every subject in a study shares
# its study effect. persons='estimated' uses plausible values -- each unit's
# mode plus a draw from its effects' posterior given its data, drawn jointly
# over the unit's levels with one seed for every level's call, so the
# differences split one draw. At a parameter draw the mode is the one
# linearised around the estimate, as everything else reporting plausible values
# at draws takes it (`ctsem_laplace_subject_values(conditional = TRUE)`); a
# joint sample that kept its effect draws pairs each parameter draw with its
# own.

# Where each level's effects sit in the raw vector: the level's declared
# `re_index`, checked at the estimate against what the engine's own subject
# vectors move. A level's own contribution to a subject's vector is the
# difference between the values built from level l outward and from level l+1
# outward, nonzero only at that level's positions.
.ctVarDecompPositions <- function(fit, levels) {
  spec <- .ctBackendSpec(fit)
  module <- .ctJuliaModule(spec$project)
  objective <- .ctJuliaObjective(fit, layer = TRUE)
  raw <- .ctJuliaNumericVector(as.numeric(fit$estimate$raw))
  nlevels <- length(levels)
  values <- lapply(seq_len(nlevels + 1L), function(l)
    as.matrix(.ctBackendJuliaValue(module$ctsem_laplace_subject_values(
      objective, raw, from_level = as.integer(l), ti_effects = FALSE))))
  lapply(seq_len(nlevels), function(l) {
    declared <- as.integer(spec$laplace$levels[[l]]$re_index)
    found <- which(apply(abs(values[[l]] - values[[l + 1L]]) > 1e-10, 2L, any))
    # An effect that is estimated at exactly zero moves nothing, so `found` can
    # be a subset of `declared`; the other way round would mean the level
    # shifts a parameter it does not declare, and this route would be placing
    # draws in the wrong coordinates.
    if (length(setdiff(found, declared))) {
      stop("The random effects of level '", levels[l], "' move raw parameters ",
        paste(setdiff(found, declared), collapse = ', '), ' that the level ',
        'does not declare, so a drawn effect cannot be placed. This is a bug ',
        '-- please report the fit.', call. = FALSE)
    }
    declared
  })
}

# Plausible values of every subject at each row of `at`, one array per level
# (draws x subjects x parameters, all levels from the same joint draws): from
# `effects`, a joint sample's saved effect draws, where given.
.ctVarDecompCoordinateDraws <- function(fit, objective, at, nlevels,
  effects = NULL) {
  module <- .ctJuliaModule(.ctBackendSpec(fit)$project)
  seed <- sample.int(.Machine$integer.max, 1L)
  nsubjects <- length(.ctBackendSpec(fit)$subject_starts)
  width <- ncol(.ctBackendRawFull(fit, at[1L, , drop = FALSE]))
  draws <- .ctJuliaPut(at)
  lapply(seq_len(nlevels), function(l) {
    value <- if (!is.null(effects)) {
      module$ctsem_laplace_subject_values(objective, draws, .ctJuliaPut(effects),
        ti_effects = FALSE, from_level = as.integer(l))
    } else module$ctsem_laplace_subject_values(objective, draws,
      .ctJuliaNumericVector(as.numeric(fit$estimate$raw)), conditional = TRUE,
      seed = as.integer(seed), from_level = as.integer(l), ti_effects = FALSE)
    array(as.numeric(.ctBackendJuliaValue(value)), c(nrow(at), nsubjects, width))
  })
}

.ctVarDecompCoordinateMats <- function(plan, fit, raw, drawrows) {
  layout <- .ctBackendSummaryLayout(fit)
  names <- .ctVarDecompLayoutNeeds(layout)
  rows <- .ctVarDecompFlatRows(layout, names)
  full <- .ctBackendRawFull(fit, raw)
  ndraws <- nrow(full)
  npeople <- length(plan$rows)
  nlevels <- plan$nlevels
  draw <- rep(seq_len(ndraws), npeople)
  person <- rep(seq_len(npeople), each = ndraws)
  population <- .ctBackendParMatricesFlat(fit, t(full))

  # Each trajectory's raw vector at each level, innermost first, built outward
  # in so that level l carries the effects of l and everything outside it.
  values <- vector('list', nlevels)
  if (identical(plan$persons, 'model')) {
    cumulative <- full[draw, , drop = FALSE]
    for (l in rev(seq_len(nlevels))) {
      at <- plan$position[[l]]
      if (length(at)) {
        covariance <- .ctBackendRawPopCovMatrix(layout, population,
          if (l == 1L) 'rawpopcov' else paste0('rawpopcov_', l))
        z <- plan$normals[[l]]
        nunits <- nrow(z)
        shift <- matrix(0, ndraws * nunits, length(at))
        for (d in seq_len(ndraws)) {
          shift[d + ndraws * (seq_len(nunits) - 1L), ] <- z %*%
            .ctSymRoot(matrix(covariance[d, , ], length(at), length(at)))
        }
        cumulative[, at] <- cumulative[, at, drop = FALSE] +
          shift[draw + ndraws * (plan$unit[person, l] - 1L), , drop = FALSE]
      }
      values[[l]] <- cumulative
    }
  } else if (nlevels) {
    objective <- .ctJuliaObjective(fit, layer = TRUE)
    replicates <- plan$replicates
    at <- raw[rep(seq_len(ndraws), replicates), , drop = FALSE]
    plausible <- .ctVarDecompCoordinateDraws(fit, objective, at, nlevels,
      if (isTRUE(plan$saved)) fit$sample$effects[drawrows, , drop = FALSE])
    pick <- draw + ndraws * (plan$copy[person] - 1L) +
      ndraws * replicates * (plan$subjects[person] - 1L)
    values <- lapply(plausible, function(p)
      matrix(p, ndraws * replicates * dim(p)[2L])[pick, , drop = FALSE])
  }

  steps <- vector('list', nlevels + 1L)
  for (s in seq_len(nlevels)) {
    out <- list()
    for (g in seq_along(plan$tipreds)) {
      at <- which(plan$group[person] == g)
      flat <- .ctBackendParMatricesFlat(fit, t(values[[s]][at, , drop = FALSE]),
        tipreds = plan$tipreds[[g]], rows = rows)
      out <- .ctVarDecompScatter(out, at, .ctVarDecompFlatArrays(flat, layout, names),
        length(draw))
    }
    steps[[s]] <- out
  }
  # The population's, at each draw and predictor row.
  out <- list()
  for (g in seq_along(plan$tipreds)) {
    flat <- if (is.null(plan$tipreds[[g]])) population[rows, , drop = FALSE] else
      .ctBackendParMatricesFlat(fit, t(full), tipreds = plan$tipreds[[g]], rows = rows)
    at <- which(plan$group[person] == g)
    out <- .ctVarDecompScatter(out, at, lapply(.ctVarDecompFlatArrays(flat, layout,
      names), function(x) x[draw[at], , , drop = FALSE]), length(draw))
  }
  steps[[nlevels + 1L]] <- out
  .ctVarDecompBatch(plan, lapply(steps, .ctVarDecompTrim, plan$nlatent,
    augmented = FALSE), ndraws)
}


# Augmented fits: a person is a carrier vector ---------------------------------
.ctVarDecompCarrierMats <- function(plan, fit, raw) {
  layout <- .ctBackendSummaryLayout(fit)
  names <- .ctVarDecompLayoutNeeds(layout)
  rows <- .ctVarDecompFlatRows(layout, names)
  full <- .ctBackendRawFull(fit, raw)
  ndraws <- nrow(full)
  npeople <- length(plan$rows)
  draw <- rep(seq_len(ndraws), npeople)
  person <- rep(seq_len(npeople), each = ndraws)
  carrier <- plan$carrier
  inner <- outer <- list()
  for (g in seq_along(plan$tipreds)) {
    at <- which(plan$group[person] == g)
    # The augmented T0 block depends on the parameters and on the predictors,
    # and on nothing else.
    population <- .ctVarDecompFlatArrays(.ctBackendParMatricesFlat(fit, t(full),
      tipreds = plan$tipreds[[g]], rows = rows), layout, names)
    outer <- .ctVarDecompScatter(outer, at, lapply(population, function(x)
      x[draw[at], , , drop = FALSE]), length(draw))
    if (!length(carrier)) next
    states <- matrix(population$T0MEANS[draw[at], , 1L], length(at))
    if (identical(plan$persons, 'model')) {
      for (d in seq_len(ndraws)) {
        here <- which(draw[at] == d)
        root <- .ctSymRoot(matrix(population$T0cov[d, carrier, carrier],
          length(carrier), length(carrier)))
        states[here, carrier] <- states[here, carrier, drop = FALSE] +
          plan$normals[person[at[here]], , drop = FALSE] %*% root
      }
    } else states[, carrier] <- plan$states[person[at], , drop = FALSE]
    flat <- .ctVarDecompAtStates(fit, t(full[draw[at], , drop = FALSE]), t(states),
      plan$tipreds[[g]], rows)
    inner <- .ctVarDecompScatter(inner, at, .ctVarDecompFlatArrays(flat, layout,
      names), length(draw))
  }
  if (!length(carrier)) inner <- outer
  inner <- .ctVarDecompTrim(inner, plan$nlatent, augmented = TRUE)
  outer <- .ctVarDecompTrim(outer, plan$nlatent, augmented = TRUE)
  .ctVarDecompBatch(plan, if (plan$nlevels) list(inner, outer) else list(inner),
    ndraws)
}

# Subject matrices at the estimate, for the carrier route's estimated persons.
.ctVarDecompSubjectMatrices <- function(fit) {
  point <- fit
  point$estimate$rawposterior <- NULL
  ctExtract(point, subjectMatrices = TRUE)
}

# A draw from each donor's stored posterior (`fit$randomeffects`), in the
# filter's state units at the carrier positions: donors x carriers. The
# posterior is stored in raw units, the state's departure over its
# state-to-raw factor, which is multiplied back here. Drawn through the
# symmetric root (`.ctSymRoot()`), which a parameter with no random effect --
# a zero row and column -- does not stop, as chol() would.
.ctVarDecompCarrierNoise <- function(fit, donors, carrier) {
  entry <- fit$randomeffects[[1L]]
  if (is.null(entry)) stop("This fit carries no posterior of its subjects' ",
    "effects (fits made before ctsem stored one); refit it, or use ",
    "persons='model'.", call. = FALSE)
  effects <- .ctBackendSpec(fit)$random_effects
  effects <- effects[effects$type %in% "sd", , drop = FALSE]
  at <- match(entry$param, as.character(effects$param))
  position <- match(as.integer(effects$row[at]), carrier)
  if (anyNA(position)) stop("The stored posterior's effects are not this fit's ",
    "carrier states. This is a bug -- please report the fit.", call. = FALSE)
  scale <- as.numeric(effects$scale[at])
  scale[!is.finite(scale) | scale == 0] <- 1
  covs <- .ctRandomEffectsCov(entry)
  out <- matrix(0, length(donors), length(carrier))
  for (index in seq_along(donors)) {
    root <- .ctSymRoot(covs[donors[index], , ])
    out[index, position] <- scale * as.numeric(root %*%
      stats::rnorm(length(entry$param)))
  }
  out
}


# Stan fits: the subject matrices the filter saved -----------------------------
#
# At each draw, the subjects' matrices the filter saved at that draw's
# parameters: estimated persons, attenuated by shrinkage. Stan allocates a
# per-subject array for a matrix exactly when that matrix's specification
# varies by subject (the `savesubjectmatrices && (sum(whenmat[..]) ||
# statedep[..])` gate in R/ctModelWriter.R), so a matrix with no `subj_` entry
# is the same for everyone and `pop_` is its value.
.ctVarDecompStanMats <- function(plan, fit, raw) {
  extracted <- suppressMessages(stan_constrainsamples(sm = fit$stanmodel,
    standata = fit$standata, samples = raw, cores = 1L, savescores = FALSE,
    savesubjectmatrices = TRUE, dokalman = TRUE, onlyfirstrow = FALSE, pcovn = 5))
  ndraws <- dim(extracted$pop_DRIFT)[1L]
  npeople <- length(plan$subjects)
  inner <- outer <- list()
  for (name in .ctVarDecompNeeded) {
    population <- extracted[[paste0('pop_', name)]]
    subject <- extracted[[paste0('subj_', name)]]
    if (is.null(population)) next
    outer[[name]] <- population[rep(seq_len(ndraws), npeople), , , drop = FALSE]
    inner[[name]] <- if (is.null(subject)) outer[[name]] else
      array(subject[, plan$subjects, , , drop = FALSE],
        c(ndraws * npeople, dim(subject)[3:4]))
  }
  inner <- .ctVarDecompTrim(inner, plan$nlatent, augmented = FALSE)
  outer <- .ctVarDecompTrim(outer, plan$nlatent, augmented = FALSE)
  .ctVarDecompBatch(plan, if (plan$nlevels) list(inner, outer) else list(inner),
    ndraws)
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
  'MANIFESTcov', 'Jy', 'THRESHOLDS', 'ASYMPTOTES', 'RESPONSE', 'WIENER')

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
# One body over whatever the person source produced, which is why the three
# representations of individual differences do not each need their own.
#
# Each trajectory's moments are collapsed over its rows, then averaged over the
# persons at each draw. The time term uses the population variance over a
# person's rows, because those rows are the design rather than a sample from
# anything; a level's between term uses the sample variance over that level's
# units, which is an estimate of a population variance.
#
# Returns, for D draws, V variables and T times (one without `times`): arrays
# D x V x T for `between`, `within.stochastic` and `within.measurement`,
# D x V for `within.deterministic`, `levels` D x levels x V x T, and the
# sources' D x V x T and D x V.
.ctVarDecompMomentComponents <- function(batch, design, levels, nlatent, latents,
  type, scale, gh, continuoustime, nmanifest, times = NULL, sources = NULL) {

  # `times` turns the window average into a decomposition at each point on a
  # common grid. Every person is then evaluated at the same times rather than
  # at their own rows, which is what makes a variance across persons at a given
  # time mean anything.
  #
  # There is no deterministic term in that form, and its absence is the point
  # rather than an omission: that term *is* the variance of the mean path over
  # t, so it exists only because the aggregate averages across the window. At
  # one instant a person's mean path is a number, not a spread.
  bytime <- !is.null(times)
  ntime <- if (bytime) length(times) else 1L
  ndraws <- batch$ndraws
  npeople <- length(batch$rows)
  nvar <- nmanifest + if (latents) nlatent else 0L
  # A grid starts at the record's start, time 0, where the initial state
  # applies; a grid that starts later steps there from 0, and that row is
  # dropped. Read as the start itself, a first time of 4 gave the initial
  # state's within person variance, which a random initial mean leaves at its
  # floor, and an ICC of 1.
  fromzero <- bytime && times[1L] > 0
  grid <- if (fromzero) c(0, times) else times
  persontimes <- lapply(batch$rows, function(rows) if (bytime) grid else
    design$time[rows])
  persontd <- if (!bytime && !is.null(design$tdpreds)) lapply(batch$rows,
    function(rows) design$tdpreds[rows, , drop = FALSE])
  response <- identical(scale, 'response') && any(type != 0L)

  value <- vector('list', length(batch$mats))
  for (s in seq_along(batch$mats)) {
    # The within person terms belong to the person's own trajectory, so they
    # come from the innermost step; the outer ones only supply the level means
    # that the between terms are differences of.
    mats <- batch$mats[[s]]
    moments <- .ctVarDecompMoments(mats, ndraws, persontimes, persontd,
      continuoustime, covariance = s == 1L || response,
      sources = if (s == 1L) sources)
    if (fromzero) moments$rows <- moments$rows[-1L]
    stats <- .ctVarDecompRowStats(moments, mats, ndraws, lengths(persontimes),
      type, scale, gh, nmanifest, nlatent, latents, bytime, within = s == 1L)
    value[[s]] <- stats$value
    if (s == 1L) inner <- stats
  }

  # Trajectories x ... to draws x ..., averaged over the persons.
  perdraw <- function(x) {
    d <- dim(x)
    colMeans(aperm(array(x, c(ndraws, npeople, d[-1L])),
      c(2L, 1L, seq_along(d)[-1L] + 1L)))
  }

  # Each level's own departure from the level outside it, varied over that
  # level's units so that a study with more subjects in it still counts once.
  contribution <- NULL
  if (length(levels)) {
    contribution <- array(0, c(ndraws, length(levels), nvar, ntime),
      dimnames = list(NULL, levels, NULL, NULL))
    for (l in seq_along(levels)) {
      difference <- matrix(aperm(array(value[[l]] - value[[l + 1L]],
        c(ndraws, npeople, nvar, ntime)), c(2L, 1L, 3L, 4L)), npeople)
      unit <- batch$unit[, l]
      means <- rowsum(difference, unit) / as.vector(rowsum(rep(1, npeople), unit))
      nunits <- nrow(means)
      contribution[, l, , ] <- if (nunits > 1L)
        colSums((means - rep(colMeans(means), each = nunits))^2) / (nunits - 1L)
      else 0
    }
  }
  out <- list(
    between = if (is.null(contribution)) array(0, c(ndraws, nvar, ntime)) else
      array(apply(contribution, c(1L, 3L, 4L), sum), c(ndraws, nvar, ntime)),
    levels = contribution,
    within.stochastic = perdraw(inner$stochastic),
    within.measurement = perdraw(inner$measurement),
    npersons = npeople)
  # The deterministic term is the variance over t and has no per-instant form.
  if (!bytime) out$within.deterministic <- perdraw(inner$deterministic)
  # Each source's share of the stochastic term, and over the window how much of
  # it is the record's mean: what does not move within a record.
  if (length(sources)) {
    partnames <- c('initial', names(sources))
    out$sources <- lapply(seq_along(partnames), function(g)
      array(perdraw(inner$parts[, , , g, drop = FALSE]), c(ndraws, nvar, ntime)))
    names(out$sources) <- partnames
    if (!bytime) {
      out$sourcelevels <- lapply(seq_along(partnames), function(g)
        matrix(perdraw(inner$partlevels[, , g, drop = FALSE]), ndraws, nvar))
      names(out$sourcelevels) <- partnames
    }
  }
  out
}

# Each trajectory's moments collapsed over its rows, through the measurement
# model: trajectories x variables x times (one without a common grid), and with
# `within` the within person terms too.
.ctVarDecompRowStats <- function(moments, mats, ndraws, rowcounts, type, scale,
  gh, nmanifest, nlatent, latents, bytime, within) {
  nvar <- nmanifest + if (latents) nlatent else 0L
  ntrajectory <- dim(mats$DRIFT)[1L]
  count <- rep(rowcounts, each = ndraws)
  ntime <- if (bytime) length(moments$rows) else 1L
  manifest <- seq_len(nmanifest)
  latent <- nmanifest + seq_len(nlatent)
  value <- array(0, c(ntrajectory, nvar, ntime))
  nparts <- if (within && !is.null(moments$levellatent))
    dim(moments$levellatent)[3L] else 0L
  if (within) {
    stochastic <- measurement <- value
    expected <- vector('list', length(moments$rows))
    if (nparts) parts <- array(0, c(ntrajectory, nvar, ntime, nparts))
  }
  for (r in seq_along(moments$rows)) {
    row <- moments$rows[[r]]
    b <- row$index
    t <- if (bytime) r else 1L
    rowvalue <- varmean <- condvar <- matrix(0, length(b), nvar)
    for (vi in manifest) {
      measured <- .ctVarDecompMeasurement(row$linearmean[, vi],
        if (!is.null(row$linearvar)) row$linearvar[, vi], type[vi],
        mats$MANIFESTcov[b, vi, vi], scale, gh)
      rowvalue[, vi] <- measured$expected
      if (within) {
        varmean[, vi] <- measured$varmean
        condvar[, vi] <- measured$condvar
      }
    }
    if (latents) {
      rowvalue[, latent] <- row$latentmean
      if (within) varmean[, latent] <- row$latentvar
    }
    value[b, , t] <- value[b, , t] + rowvalue
    if (within) {
      stochastic[b, , t] <- stochastic[b, , t] + varmean
      measurement[b, , t] <- measurement[b, , t] + condvar
      expected[[r]] <- rowvalue
      if (nparts) {
        parts[b, manifest, t, ] <- parts[b, manifest, t, ] + as.vector(row$partlinear)
        if (latents) parts[b, latent, t, ] <- parts[b, latent, t, ] +
          as.vector(row$partlatent)
      }
    }
  }
  out <- list()
  if (!bytime) {
    value[, , 1L] <- value[, , 1L] / count
    if (within) {
      stochastic[, , 1L] <- stochastic[, , 1L] / count
      measurement[, , 1L] <- measurement[, , 1L] / count
      if (nparts) for (g in seq_len(nparts)) parts[, , 1L, g] <- parts[, , 1L, g] / count
      # The variance of the expected value over the person's own rows.
      mean <- matrix(value[, , 1L], ntrajectory)
      deterministic <- matrix(0, ntrajectory, nvar)
      for (r in seq_along(moments$rows)) {
        b <- moments$rows[[r]]$index
        deterministic[b, ] <- deterministic[b, , drop = FALSE] +
          (expected[[r]] - mean[b, , drop = FALSE])^2
      }
      out$deterministic <- deterministic / count
    }
  }
  out$value <- value
  if (within) {
    out$stochastic <- stochastic
    out$measurement <- measurement
    if (nparts) {
      out$parts <- parts
      out$partlevels <- array(0, c(ntrajectory, nvar, nparts))
      out$partlevels[, manifest, ] <- moments$levellinear
      if (latents) out$partlevels[, latent, ] <- moments$levellatent
    }
  }
  out
}

# Chunks' components, bound along their draws.
.ctVarDecompBindDraws <- function(chunks) {
  if (length(chunks) == 1L) return(chunks[[1L]])
  bind <- function(parts) {
    first <- parts[[1L]]
    if (is.null(first)) return(NULL)
    if (is.list(first)) {
      out <- lapply(seq_along(first), function(i) bind(lapply(parts, `[[`, i)))
      names(out) <- names(first)
      return(out)
    }
    if (is.null(dim(first))) return(first)
    stacked <- do.call(rbind, lapply(parts, function(x) matrix(x, dim(x)[1L])))
    array(stacked, c(nrow(stacked), dim(first)[-1L]),
      dimnames = if (!is.null(dimnames(first))) c(list(NULL), dimnames(first)[-1L]))
  }
  bind(chunks)
}


# The simulation route --------------------------------------------------------
#
# For a model whose dynamics depend on the state there is no moment recursion to
# run: the transition is a different transition at every point the trajectory
# visits. The engine draws the trajectory instead, through the same state pass
# `ctGenerate(intoverstates = FALSE)` uses, and the decomposition is taken over
# the draws. It is computed at the estimate only: at every parameter draw it
# would be a simulation of its own.
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
          linear <- matrix(mats$MANIFESTMEANS[1L, vi, 1L], npaths, length(rows))
          for (li in seq_len(nlatent)) {
            linear <- linear + mats$LAMBDA[1L, vi, li] * dynamic[[li]]
          }
          linear
        }
        conditional <- if (latentonly) 0 else .ctVarDecompSimulatedCondVar(g,
          type[vi], mats$MANIFESTcov[1L, vi, vi], scale)
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

  # As the moment route returns them, for the one draw that is the estimate.
  single <- function(x) array(x, c(1L, nvar, 1L))
  list(
    between = single(if (npeople > 1L)
      pmax(apply(personmean, 2L, stats::var) - colMeans(betweencorrection), 0)
      else rep(0, nvar)),
    within.deterministic = matrix(pmax(colMeans(persondet) -
      colMeans(detcorrection), 0), 1L),
    within.stochastic = single(colMeans(personstoch)),
    within.measurement = single(colMeans(personmeas)),
    npersons = npeople)
}

# One person's matrices, materialised at a whole filter state, as arrays with
# a single leading entry. julia only. Messages suppressed, and said once by the
# print method instead: called once per person, each call otherwise announces
# which cells are state dependent and where they were evaluated.
.ctVarDecompPersonAt <- function(fit, tipreds, state, nlatent) {
  # filterstate=, not state=: the carrier entries are this person's random
  # effects and are the whole point of materialising here.
  drawn <- suppressMessages(ctBackendParMatrices(fit, tipreds = tipreds,
    filterstate = state, trim = FALSE))
  mats <- list()
  for (name in .ctVarDecompNeeded) {
    if (is.null(drawn[[name]])) next
    value <- as.matrix(drawn[[name]])
    mats[[name]] <- array(value, c(1L, dim(value)))
  }
  .ctVarDecompTrim(mats, nlatent, augmented = TRUE)
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


# The returned table ----------------------------------------------------------
#
# Every column at every draw, draws x rows, from the components: the parts,
# then what is built from them, each built at its own draw.
.ctVarDecompColumns <- function(components, parts) {
  flat <- function(x) matrix(x, dim(x)[1L])
  columns <- lapply(stats::setNames(parts, parts), function(part)
    flat(components[[part]]))
  # One column per random effect level, where there is more than one: with a
  # single level it would repeat `between` under another name.
  levels <- components$levels
  if (!is.null(levels) && dim(levels)[2L] > 1L) {
    for (l in seq_len(dim(levels)[2L])) {
      columns[[paste0('between.', dimnames(levels)[[2L]][l])]] <-
        flat(levels[, l, , , drop = FALSE])
    }
  }
  # The stochastic term by source, and the part of each that is the record's
  # mean, after the parts they split.
  for (g in names(components$sources)) {
    columns[[paste0('within.stochastic.', g)]] <- flat(components$sources[[g]])
  }
  for (g in names(components$sourcelevels)) {
    columns[[paste0('level.', g)]] <- flat(components$sourcelevels[[g]])
  }
  columns$within <- Reduce(`+`, columns[setdiff(parts, 'between')])
  columns$total <- columns$between + columns$within
  # Each part followed by the parts it splits: the levels after `between`, the
  # sources and their record levels after `within.stochastic`.
  levelcols <- grep('^between\\.', names(columns), value = TRUE)
  variances <- c('between', levelcols,
    setdiff(parts, c('between', 'within.measurement')),
    grep('^(within\\.stochastic\\.|level\\.)', names(columns), value = TRUE),
    'within.measurement', 'within', 'total')
  # A variable with no variance at all has no proportions, and NA says that
  # where 0/0 would print NaN. A constant auxiliary state is the case: a
  # higher order model's carried coordinate contributes nothing to anything.
  share <- function(part, whole) ifelse(whole > 0, part / whole, NA_real_)
  for (part in setdiff(variances, 'total')) {
    columns[[paste0('prop.', part)]] <- share(columns[[part]], columns$total)
  }
  # Intraclass correlations. `icc` is the between person share of the total,
  # and `icc.process` the same with measurement error left out of the total. A
  # level's is the share of that level and every level above it: what two
  # observations sharing a unit of that level and no inner one have in common
  # (Snijders & Bosker's definition), so the innermost level's is `icc`.
  columns$icc <- share(columns$between, columns$total)
  columns$icc.process <- share(columns$between,
    columns$total - columns$within.measurement)
  if ('between.process' %in% levelcols) stop("A random effect level is called ",
    "'process', which would name its intraclass correlation as icc.process does. ",
    'Rename that id column.', call. = FALSE)
  for (l in seq_along(levelcols)) {
    columns[[sub('^between', 'icc', levelcols[l])]] <- share(
      Reduce(`+`, columns[levelcols[l:length(levelcols)]]), columns$total)
  }
  order <- c(variances, paste0('prop.', setdiff(variances, 'total')), 'icc',
    'icc.process', sub('^between', 'icc', levelcols))
  columns[order]
}


#' Decompose model implied variance into between person, process and measurement
#'
#' Splits the variance each indicator of a fitted model implies into a between
#' person part, a within person part from the latent process, and a within
#' person part from measurement error. The process part is further split into
#' the contribution of the system noise and that of the mean path changing
#' over time. Each part is also given as a proportion of the total, along with
#' intraclass correlations. Each value is computed at a set of the fit's
#' parameter draws, giving intervals.
#'
#' @param fit fit object as generated by \code{\link{ctFit}}, from either
#'   backend.
#' @param method How the model implied moments are computed. \code{'moment'}
#'   uses the forward moment recursion, which is exact and needs linear
#'   dynamics. \code{'simulation'} draws trajectories with the julia engine,
#'   for models with state dependent \code{DRIFT} or \code{DIFFUSION} cells.
#'   \code{'auto'}, the default, uses \code{'moment'} unless the model has
#'   such cells.
#' @param npaths Trajectories drawn per person when \code{method='simulation'}.
#'   The result is unbiased for any \code{npaths}, but small values give noisy
#'   results.
#' @param persons Which persons the between person variance refers to.
#'   \code{'model'} draws persons from the fitted population distribution of
#'   the random effects. \code{'estimated'} uses the subjects in the data. On a
#'   \code{backend='julia'} fit each subject is represented by its estimated
#'   effects plus a draw from their posterior. On a stan fit each subject is
#'   represented by its estimated (empirical Bayes) matrices, which understates
#'   the between person variance. \code{'auto'}, the default, uses
#'   \code{'model'} for julia fits and \code{'estimated'} for stan fits, and
#'   says which it used.
#' @param scale For a binary or ordinal indicator, whether to decompose the
#'   linear predictor behind the logit (\code{'latent'}, the default, with a
#'   measurement variance of \eqn{\pi^2/3}) or the observed response
#'   (\code{'response'}, integrating the link over each row's normal
#'   distribution). The two are the same for a Gaussian indicator.
#' @param times If given, decompose at each of these times. By default the
#'   decomposition averages over each person's observation times. Times are
#'   measured from the start of each person's record, where the initial state
#'   applies. Either a numeric vector, or \code{'asdata'} for the times of
#'   the observed occasions since each record began. Every person is
#'   evaluated at the same times. The result gains a \code{time} column and
#'   has no \code{within.deterministic} column, since that term is the variance
#'   of the mean path over time. Not available with \code{method='simulation'}
#'   or with time dependent predictors.
#' @param npersons Number of persons for the between person variance, at each
#'   parameter draw. Each subject's observation times are used
#'   \code{round(npersons / N)} times for N subjects, and at least once. With
#'   \code{persons='model'} each use gets its own draw of random effects from
#'   the population distribution, and with \code{persons='estimated'} on a
#'   julia fit its own draw from the subject's posterior.
#' @param nsamples Number of parameter draws to compute the decomposition at.
#'   Each value is the median of that quantity over the draws, and its 2.5\%
#'   and 97.5\% quantiles give an interval. 0 computes the decomposition once
#'   at the estimate, without intervals. \code{method='simulation'}, and
#'   \code{persons='estimated'} on an \code{intoverpop='augmented'} fit, give
#'   values at the estimate only.
#' @param latents If TRUE, also decompose the latent processes, which have no
#'   measurement component.
#' @param quadpoints Gauss-Hermite nodes used when \code{scale='response'}.
#' @param subjects \code{'all'}, or the subjects whose designs the decomposition
#'   is computed over, as ids in the data (see \code{realid}).
#' @param realid If TRUE (the default), subjects are read as the ids in the data,
#' and integers that are not ids there are read as ctsem's numbering of the
#' subjects, 1 to N, with a message. If FALSE, they are read as that numbering.
#' @param sources Optional named list of latent names, one entry per group of
#'   processes whose innovations are independent of the other groups, such as
#'   fast affect, a daily process and a slow process. The stochastic within
#'   person term is then split by source, into \code{within.stochastic.<group>}
#'   for each group's innovations and \code{within.stochastic.initial} for
#'   what remains of the initial state variance. These add up to
#'   \code{within.stochastic}. When averaging over the observation times,
#'   \code{level.<group>} gives the part of each group's term that is variance
#'   in a person's average level over their record. This shows how much a slow
#'   process acts as a stable level within the record. Needs
#'   \code{method='moment'}.
#'
#' @details The decomposition applies the law of total variance over the
#'   (person, occasion) pairs in the design of the fit. With
#'   \eqn{g = E[y | \eta]} the expected observation given the latent state,
#'
#'   \deqn{Var(y) = E[Var(y|\eta)] + E_i E_t Var(g) + E_i Var_t(E[g]) + Var_i(E[g])}
#'
#'   and the four terms are the four columns reported. They are non-negative
#'   and sum to the total.
#'
#'   \code{within.deterministic} is the variance of a person's mean path over
#'   their observation times. It is zero for a stationary process without time
#'   dependent predictors, and positive when there is a trend, a shift due to
#'   a predictor, or a process that has not settled from its initial
#'   distribution.
#'
#'   The moments are computed over each subject's observation times without
#'   using the observed values, so unequal spacing, missing observations and
#'   time dependent predictors are all accounted for. The result is exact for
#'   models with linear dynamics. For a model whose \code{DRIFT} or
#'   \code{DIFFUSION} depends on the state, \code{method='moment'} is refused,
#'   and \code{method='simulation'} applies the same decomposition to simulated
#'   trajectories. State dependent measurement cells (\code{LAMBDA},
#'   \code{MANIFESTMEANS}, \code{MANIFESTVAR}) are not supported by either
#'   method.
#'
#'   \strong{Intervals.} The decomposition is computed at each of
#'   \code{nsamples} draws of the fit's parameters: draws from the normal
#'   approximation for an optimised fit, and posterior draws for a sampled
#'   one. Each value in the returned table is the median of that quantity over
#'   the draws. Each median is taken separately, so the medians of the parts
#'   need not add up to the median of the total, though at each draw they add
#'   up exactly. The same persons are used at every draw.
#'
#'   \strong{Random effect levels.} A model with a grouping level above the
#'   subject (\code{id = c('subject','study')} in \code{\link{ctModel}}),
#'   fitted with \code{intoverpop='laplace'} or sampled, gets one
#'   \code{between.<level>} column per level. Each is the variance of that
#'   level's departure from the level above it, and together they sum to
#'   \code{between}. Estimated random effects are shrunk towards the level
#'   above, more strongly where a unit has little data, so
#'   \code{persons='estimated'} adds to each estimate a draw from that unit's
#'   posterior.
#'
#'   \strong{When the split changes over time.} The columns are averages over
#'   each person's occasions. For a process that has not settled the split
#'   differs between occasions, and \code{within.deterministic} is then
#'   large. Use \code{times} to see the split at each time.
#'
#'   \strong{Designs.} Each person is evaluated at their own observation
#'   times, so balanced and unbalanced designs and missing observations are
#'   all handled, and the result does not depend on which observations are
#'   missing. Persons are weighted equally, however many observations each
#'   has.
#'
#'   With \code{method='simulation'}, the function checks that the engine
#'   draws the individually varying parameters, and stops with an error if it
#'   does not.
#'
#'   \strong{Intraclass correlations.} \code{icc} is the between person
#'   variance divided by the total, and \code{icc.process} the between person
#'   variance divided by the total without measurement error. For a latent
#'   process the two are the same. With more than one random effect level,
#'   \code{icc.<level>} is the between person variance of that level and all
#'   levels above it, divided by the total. For a level above the person, such
#'   as a study, this is the correlation between observations of two persons
#'   in the same study. For the innermost level it equals \code{icc}.
#'
#' @return A data frame of class \code{ctVarianceDecomposition}, one row per
#'   variable, with columns \code{variable}, \code{type}, \code{between},
#'   \code{within.deterministic}, \code{within.stochastic},
#'   \code{within.measurement}, \code{within} and \code{total}, the
#'   \code{between.<level>}, \code{within.stochastic.<source>} and
#'   \code{level.<source>} columns where they apply, a \code{prop.} column for
#'   each of these giving it as a proportion of \code{total}, and the
#'   intraclass correlations \code{icc}, \code{icc.process} and
#'   \code{icc.<level>}. With \code{nsamples} above 0 each value is a median
#'   over parameter draws, \code{attr(x, 'intervals')} is a data frame with
#'   the 2.5\%, 50\% and 97.5\% quantiles of every column for every variable,
#'   and \code{attr(x, 'draws')} holds the values at each draw, an array of
#'   draws by rows by columns. The population the between column refers to,
#'   the scale, and the number of persons and draws are recorded as attributes
#'   and shown by \code{print}.
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
  times = NULL, npersons = 200L, npaths = 20L, nsamples = 1000L, latents = TRUE,
  quadpoints = 21L, subjects = 'all', realid = TRUE, sources = NULL) {

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
  if (!is.null(sources)) {
    if (!identical(method, 'moment')) stop("sources needs method='moment': the ",
      'split follows the covariance recursion, which the simulation route does ',
      'not run.', call. = FALSE)
    if (identical(scale, 'response') &&
        any(as.integer(.ctFitModelObject(fit)$manifesttype) != 0L)) stop(
      "sources splits variance on the latent scale, where it adds up; use ",
      "scale='latent' with non-Gaussian indicators.", call. = FALSE)
  }
  if (identical(method, 'simulation') && !.ctFitIsJulia(fit)) {
    stop("method='simulation' needs a backend='julia' fit: it draws the ",
      'trajectories from the engine, which the stan path has no entry point ',
      'for. Refit with backend=\'julia\'.', call. = FALSE)
  }
  if (identical(method, 'simulation') && .ctFitIsJulia(fit) &&
      .ctSpecEffectsAreCoordinates(.ctBackendSpec(fit))) {
    stop("method='simulation' is not available when the random effects are ",
      'separate coordinates: they are not carrier states, so there is no part ',
      'of the innovation draw that pins a person while its path is redrawn.',
      call. = FALSE)
  }

  model <- .ctFitModelObject(fit)
  nlatent <- .ctFitNlatent(fit)
  continuoustime <- isTRUE(model$continuoustime)
  manifestNames <- model$manifestNames
  latentNames <- model$latentNames[seq_len(nlatent)]
  type <- as.integer(model$manifesttype)

  # Whether the random effects are separate coordinates, which is true of
  # `intoverpop='laplace'` and of `intoverpop='none'` -- the route a sampled fit
  # with random effects takes. Both need the level route below; they differ in
  # whether a subject's effect is a mode or a draw, which the message names.
  coordinates <- .ctFitIsJulia(fit) &&
    .ctSpecEffectsAreCoordinates(.ctBackendSpec(fit))
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
  if (!is.null(times)) {
    if (identical(method, 'simulation')) {
      stop("times= is not available for method='simulation': the drawn ",
        'trajectories are over the observed rows, so there is no common grid ',
        'to evaluate every person on.', call. = FALSE)
    }
    if (!is.null(design$tdpreds)) {
      stop('times= is not available for a model with time dependent ',
        'predictors: a grid of times carries no predictor values, and ',
        'evaluating one person at another one\'s would not be that person. ',
        'Drop times= for the window average over the observed design.',
        call. = FALSE)
    }
    # Times run from the start of each record, where the model's initial
    # state applies; 'asdata' gives each observation's time since its own
    # record began.
    times <- if (identical(as.character(times)[1L], 'asdata')) {
      start <- tapply(design$time, design$subject, min)
      sort(unique(round(design$time - start[as.character(design$subject)], 10)))
    } else sort(as.numeric(times))
    if (length(times) < 1L) stop('times= is empty.', call. = FALSE)
    if (any(!is.finite(times) | times < 0)) stop('times= are measured from ',
      'the start of each record, so must be 0 or more.', call. = FALSE)
  }
  wanted <- if (identical(subjects, 'all')) seq_len(design$nsubjects) else
    as.vector(.ctResolveSubjects(fit, subjects, realid))
  wanted <- wanted[wanted %in% unique(design$subject)]
  if (!length(wanted)) stop('No rows for the requested subjects.', call. = FALSE)

  nmanifest <- length(manifestNames)
  levels <- vapply(.ctFitRandomEffectLevels(fit), function(x) x$name,
    character(1L))

  # The parameter draws. Two routes compute at the estimate only, and say so.
  carrierestimated <- .ctFitIsJulia(fit) && !coordinates &&
    identical(persons, 'estimated') &&
    length(.ctVarDecompCarrier(fit, nlatent)) > 0L
  draws <- .ctVarDecompDraws(fit, nsamples, pointonly =
    if (identical(method, 'simulation')) paste0("Values are at the estimate, ",
      "without intervals: method='simulation' would need a simulation at every ",
      'parameter draw.')
    else if (carrierestimated) paste0("Values are at the estimate, without ",
      "intervals: persons='estimated' on an augmented fit draws each subject's ",
      'effects from a posterior the fit stored at its estimate.'))

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
  } else {
    plan <- .ctVarDecompPlan(fit, design, nlatent, persons, as.integer(npersons),
      wanted, levels, draws$rows)
    gh <- if (identical(scale, 'response')) .ctVarDecompGaussHermite(quadpoints) else NULL
    # Draws in chunks small enough that every person's matrices at every draw
    # of a chunk stay a few tens of megabytes.
    width <- max(nlatent, nmanifest, length(plan$carrier) + nlatent)
    per <- max(1L, floor(2e6 / (length(plan$rows) * width^2 *
      (length(levels) + 1L))))
    ndraws <- nrow(draws$raw)
    chunks <- split(seq_len(ndraws), ceiling(seq_len(ndraws) / per))
    sources <- .ctVarDecompSources(fit, sources, nlatent)
    .ctVarDecompBindDraws(lapply(chunks, function(index) {
      batch <- .ctVarDecompMaterialise(plan, fit, draws$raw[index, , drop = FALSE],
        draws$rows[index])
      .ctVarDecompMomentComponents(batch, design, levels, nlatent, latents, type,
        scale, gh, continuoustime, nmanifest, times = times, sources = sources)
    }))
  }

  variables <- c(manifestNames, if (latents) latentNames)
  kinds <- c(rep('manifest', nmanifest), if (latents) rep('latent', nlatent))
  bytime <- !is.null(times)
  parts <- if (bytime) c('between', 'within.stochastic', 'within.measurement')
    else c('between', 'within.deterministic', 'within.stochastic',
      'within.measurement')
  columns <- .ctVarDecompColumns(components, parts)

  # A draw the model cannot be evaluated at -- an explosive drift far in a
  # normal approximation's tail, say -- is left out of the draws, and said.
  if (isTRUE(draws$intervals)) {
    usable <- apply(is.finite(columns$total), 1L, all)
    if (!any(usable)) stop('No parameter draw gave a decomposition. Use ',
      'nsamples = 0 for the estimate alone.', call. = FALSE)
    if (!all(usable)) {
      message(sum(!usable), ' of ', length(usable), ' parameter draws gave no ',
        'decomposition and were left out.')
      columns <- lapply(columns, function(x) x[usable, , drop = FALSE])
    }
  }

  out <- if (!bytime) {
    data.frame(variable = variables, type = kinds, stringsAsFactors = FALSE)
  } else {
    # One row per (time, variable), times varying slowest so that each block of
    # the printed frame is one instant. The components run variables fastest
    # within a time, so column major order is already the row order --
    # transposing it interleaves the two and produces a table where every row
    # sums correctly and no row is about the variable it names.
    data.frame(time = rep(times, each = length(variables)),
      variable = rep(variables, length(times)),
      type = rep(kinds, length(times)), stringsAsFactors = FALSE)
  }
  keys <- out
  for (name in names(columns)) {
    value <- columns[[name]]
    out[[name]] <- if (nrow(value) == 1L) as.numeric(value) else
      apply(value, 2L, stats::median)
  }
  rownames(out) <- NULL

  if (isTRUE(draws$intervals)) {
    probs <- c(0.025, 0.5, 0.975)
    intervals <- do.call(rbind, lapply(names(columns), function(name) {
      quantiles <- matrix(t(apply(columns[[name]], 2L, stats::quantile,
        probs = probs, na.rm = TRUE, names = FALSE)), ncol = length(probs))
      colnames(quantiles) <- paste0(100 * probs, '%')
      cbind(keys, quantity = name, as.data.frame(quantiles, check.names = FALSE),
        stringsAsFactors = FALSE)
    }))
    rownames(intervals) <- NULL
    attr(out, 'intervals') <- intervals
    attr(out, 'draws') <- array(unlist(columns, use.names = FALSE),
      c(nrow(columns$total), nrow(out), length(columns)),
      dimnames = list(NULL, NULL, names(columns)))
  }
  attr(out, 'nsamples') <- if (isTRUE(draws$intervals)) nrow(columns$total) else 0L
  attr(out, 'method') <- method
  attr(out, 'persons') <- persons
  # The julia routes give estimated persons at plausible values; stan's are its
  # modes.
  attr(out, 'plausible') <- identical(persons, 'estimated') && .ctFitIsJulia(fit)
  attr(out, 'scale') <- scale
  attr(out, 'nongaussian') <- any(type != 0L)
  attr(out, 'npersons') <- components$npersons
  attr(out, 'npaths') <- if (identical(method, 'simulation')) as.integer(npaths)
  attr(out, 'continuoustime') <- continuoustime
  attr(out, 'bytime') <- bytime
  class(out) <- c('ctVarianceDecomposition', 'data.frame')
  out
}


.ctVarDecompPersonsLabel <- function(x) {
  n <- attr(x, 'npersons')
  if (identical(attr(x, 'persons'), 'model'))
    paste0(n, ' persons drawn from the fitted population distribution')
  else if (isTRUE(attr(x, 'plausible')))
    paste0(n, ' plausible values of the subjects in the data')
  else paste0('the ', n, ' estimated subjects, attenuated by shrinkage')
}

#' @export
print.ctVarianceDecomposition <- function(x, digits = 3L, ...) {
  bytime <- isTRUE(attr(x, 'bytime'))
  frame <- as.data.frame(x)
  # The named columns of `table`, rounded, beside the time, variable and type
  # each row is about: a manifest and a latent of one name are two rows.
  show <- function(table, columns, labels = columns) {
    keys <- intersect(c('time', 'variable', 'type'), names(table))
    table <- table[c(keys, columns)]
    numeric <- vapply(table, is.numeric, logical(1L))
    table[numeric] <- lapply(table[numeric], round, digits)
    names(table) <- c(keys, labels)
    print(table, row.names = FALSE)
  }
  parts <- intersect(c('between', 'within.deterministic', 'within.stochastic',
    'within.measurement'), names(frame))
  levelcols <- grep('^between\\.', names(frame), value = TRUE)
  sourcecols <- grep('^within\\.stochastic\\.', names(frame), value = TRUE)
  recordcols <- grep('^level\\.', names(frame), value = TRUE)
  # icc.process is icc where there is no measurement error at all.
  measured <- any(frame$within.measurement > 0, na.rm = TRUE)
  iccs <- c('icc', if (measured) 'icc.process', sub('^between', 'icc', levelcols))
  intervals <- attr(x, 'intervals')
  nsamples <- attr(x, 'nsamples')

  cat('Model implied variance decomposition\n\n')
  show(frame, c(parts, 'total'))
  cat('\nProportions of the total\n')
  show(frame, paste0('prop.', c(parts, 'within')), c(parts, 'within'))
  if (length(levelcols)) {
    cat('\nBetween person variance by level, and as a proportion of the total\n')
    show(frame, c(levelcols, paste0('prop.', levelcols)))
  }
  if (length(sourcecols)) {
    cat('\nStochastic within person variance by source, as proportions of the total\n')
    show(frame, paste0('prop.', sourcecols), sub('^within\\.stochastic\\.', '', sourcecols))
  }
  if (length(recordcols)) {
    cat("\nOf which in a person's average over their record\n")
    show(frame, paste0('prop.', recordcols), sub('^level\\.', '', recordcols))
  }
  if (is.null(intervals)) {
    cat('\nIntraclass correlations\n')
    show(frame, iccs)
  } else {
    cat('\nIntraclass correlations, medians and 95% intervals\n')
    for (icc in iccs) {
      cat('\n', icc, '\n', sep = '')
      show(intervals[intervals$quantity == icc, , drop = FALSE],
        c('2.5%', '50%', '97.5%'))
    }
  }

  cat('\nicc is the proportion of the total between persons',
    if (measured) ', and icc.process the same proportion with measurement error left out of the total',
    '.\n', sep = '')
  if (length(levelcols)) cat(paste(iccs[-seq_len(1L + measured)], collapse = ', '),
    ' include the between person variance of their own level and every level ',
    'above it.\n', sep = '')
  if (isTRUE(nsamples > 0L)) {
    cat('Values are medians over ', nsamples, ' parameter draws, each ',
      'quantity taken separately, so the parts need not add up to the total ',
      'exactly.\n', sep = '')
  } else cat('Values are at the estimate, without intervals.\n')
  if (bytime) {
    cat('At one time point there is no deterministic within person term: it ',
      'is the variance of the mean path over time, so it exists only in the ',
      'window average. Call without times= for that.\n', sep = '')
    cat('Between person variance refers to ', .ctVarDecompPersonsLabel(x),
      '.\n', sep = '')
    return(invisible(x))
  }
  cat('Computed by the ', attr(x, 'method'), ' route',
    if (identical(attr(x, 'method'), 'simulation'))
      paste0(', ', attr(x, 'npaths'), ' paths per person') else '', '.\n', sep = '')
  cat('Between person variance refers to ', .ctVarDecompPersonsLabel(x),
    if (isTRUE(nsamples > 0L)) ' at each draw', '.\n', sep = '')
  cat('Within person variance is over each person\'s own observation times; ',
    'the deterministic part is the mean path moving.\n', sep = '')
  if (identical(attr(x, 'method'), 'simulation')) {
    cat('The dynamics are integrated along drawn trajectories rather than ',
      'evaluated at one state, so no evaluation point is reported; the ',
      'measurement model does not depend on the state, which was checked.\n',
      sep = '')
  }
  if (isTRUE(attr(x, 'nongaussian'))) {
    if (identical(attr(x, 'scale'), 'latent')) {
      cat("Non-Gaussian indicators are decomposed on the latent response scale ",
        '(measurement variance pi^2/3).\n', sep = '')
    } else cat('Non-Gaussian indicators are decomposed on the response scale.\n')
  }
  invisible(x)
}
