# A coreset for importance sampling on the Laplace route.
#
# A Laplace fit's log density is a sum of per-unit terms, each an inner mode
# solve, and an importance sampler pays for every unit at every draw. Over the
# region the sampler visits, the units' terms vary in about as many dimensions
# as the sampler moves, plus curvature -- to first order a unit's term is its
# score times the displacement -- however many units there are. So M of the N
# units with non-negative weights and a free intercept, fitted by least
# squares to the full sum, follow it, and the sampler can pay M/N of a full
# evaluation per draw. On the strong-shell model of
# review/UNCERTAINTY-coreset-test-2026-10-21.md, 30-40 units carried the sum of
# 100, 250 or 1,000 along 'checked''s directions, and the sampler reached the
# full sum's posterior in 2.7, 4.9 and 10 times fewer evaluations.
#
# `.ctCoresetLpg()` wraps the sampler's density, and every decision it makes is
# by the cost of a sampler round in full evaluations. It evaluates rounds in
# full, as the sampler would without it, and learns from them until some
# coreset is expected to cost less than a full round; then it evaluates that
# coreset's units only (`ctsem_coreset_batch`, on the pool a full evaluation
# uses), checks a few of every round's draws in full, and evaluates a round in
# full when its check fails. The size is the one whose expected round cost is
# least (`.ctCoresetChoose()`): its share of the round's values and gradients,
# its checks, and the chance a check fails times the full round that costs.
# That chance is read from held-out errors and then from the failures the
# coreset actually has; when they make it dearer than a full round it goes
# back to full rounds and tries again with what they teach. At the end it
# checks the draws that carry the weight and corrects them where needed, so
# what it hands back is the full density's.

# The tolerances, in log density units. A weight is off by exp(error), so an
# RMS of 0.1 costs about 1% of an effective sample; `max` keeps any one draw
# within a factor of 1.6. A round's checks may drift further before the
# coreset is rebuilt (`drift`), so that it is not rebuilt on noise.
#
# Both are read where a draw can carry weight. A sampler's first rounds are
# mostly tails, 50-150 log units below its best draw, where a coreset that is
# within 0.1 near the top errs by a unit or two -- and a draw 50 down carries
# exp(-50) of the best one's weight. So a draw counts fully within `band` of
# the best and e-fold less every `scale` below that (`.ctCoresetRelevance()`),
# in the fit as in the test, and `max` reads only draws of relevance at least
# `relevant`, about 55 below the best.
.ctCoresetTol <- list(rms = 0.1, max = 0.5, drift = 1.5, final = 0.5,
  band = 10, scale = 20, relevant = 0.1)

.ctCoresetRelevance <- function(y, top = max(y)) {
  exp(-pmax(0, top - y - .ctCoresetTol$band) / .ctCoresetTol$scale)
}

# Whether errors `err` at draws of relevance `rel` pass the tolerances, `mult`
# times over.
.ctCoresetPasses <- function(err, rel, mult = 1) {
  if (!length(err) || sum(rel) <= 0) return(TRUE)
  rms <- sqrt(sum(rel * err^2) / sum(rel))
  near <- rel >= .ctCoresetTol$relevant
  worst <- if (any(near)) max(abs(err[near])) else 0
  rms <= mult * .ctCoresetTol$rms && worst <= mult * .ctCoresetTol$max
}

# What a sampler round on a coreset of `M` of `N` units is expected to cost, in
# full evaluations: its share of the round's `values` and `grads` (a gradient
# as two), the `checks` evaluated in full, and -- with probability `fail` --
# the rest of the round evaluated in full after a failed check. A full round
# costs `values + 2 * grads`.
.ctCoresetRoundCost <- function(M, N, round, fail) {
  M / N * (round$values + 2 * round$grads) + round$checks +
    fail * (round$values - round$checks)
}

# How often a round's check would fail, from errors `err` at draws of
# relevance `rel`: rounds of `checks` draws resampled in proportion to
# relevance, as a round's checks lean towards where the weight sits. On a
# private random stream, so that the sampler's own is the same with a coreset
# as without.
.ctCoresetFailRate <- function(err, rel, checks) {
  if (!length(err) || sum(rel) <= 0) return(0)
  had <- exists('.Random.seed', envir = globalenv(), inherits = FALSE)
  if (had) old <- get('.Random.seed', envir = globalenv(), inherits = FALSE)
  on.exit(if (had) assign('.Random.seed', old, envir = globalenv()) else
    if (exists('.Random.seed', envir = globalenv(), inherits = FALSE))
      rm('.Random.seed', envir = globalenv()))
  set.seed(20261021L)
  p <- rel / sum(rel)
  n <- min(checks, length(err))
  mean(vapply(seq_len(200L), function(b) {
    i <- sample.int(length(err), n, replace = TRUE, prob = p)
    !.ctCoresetPasses(err[i], rel[i], .ctCoresetTol$drift)
  }, logical(1)))
}

# Lawson-Hanson non-negative least squares.
.ctCoresetNnls <- function(A, b) {
  n <- ncol(A)
  x <- numeric(n)
  P <- logical(n)
  w <- drop(crossprod(A, b))
  tol <- 1e-10 * max(1, max(abs(w)))
  for (it in seq_len(3 * n)) {
    if (!any(!P & w > tol)) break
    P[which.max(ifelse(P, -Inf, w))] <- TRUE
    repeat {
      z <- numeric(n)
      z[P] <- qr.coef(qr(A[, P, drop = FALSE]), b)
      z[is.na(z)] <- 0
      if (all(z[P] > 0)) break
      neg <- P & z <= 0
      alpha <- min(x[neg] / (x[neg] - z[neg]))
      x <- x + alpha * (z - x)
      P <- P & x > 1e-12
      x[!P] <- 0
    }
    x <- z
    w <- drop(crossprod(A, b - A %*% x))
  }
  x
}

# Weights on units `S` for the training `rows` of `X` (draws x units) and `y`
# (the full sum), the intercept free: least squares on the centred rows, each
# weighted by its relevance (`.ctCoresetRelevance()`).
.ctCoresetWeights <- function(X, y, S, rows, rel = .ctCoresetRelevance(y[rows])) {
  Xs <- X[rows, S, drop = FALSE]
  share <- rel / sum(rel)
  mx <- colSums(Xs * share)
  my <- sum(y[rows] * share)
  root <- sqrt(rel)
  w <- .ctCoresetNnls(sweep(Xs, 2, mx) * root, (y[rows] - my) * root)
  list(S = S, w = w, c = my - sum(mx * w))
}
.ctCoresetPredict <- function(cf, Xrows) drop(cf$c + Xrows[, cf$S, drop = FALSE] %*% cf$w)

# Forward selection to `M` units from `start`: each step adds the unit whose
# term, orthogonalised against the units already in, best explains what is
# left of the centred sum, and refits the weights. Nested, so a path serves
# every size as a prefix.
.ctCoresetPath <- function(X, y, M, rows, start = integer(),
  rel = .ctCoresetRelevance(y[rows])) {
  share <- rel / sum(rel)
  root <- sqrt(rel)
  Xr <- X[rows, , drop = FALSE]
  Xc <- sweep(Xr, 2, colSums(Xr * share)) * root
  yc <- (y[rows] - sum(y[rows] * share)) * root
  S <- start
  r <- if (length(S)) yc - drop(Xc[, S, drop = FALSE] %*%
    .ctCoresetNnls(Xc[, S, drop = FALSE], yc)) else yc
  while (length(S) < min(M, ncol(X))) {
    cand <- setdiff(seq_len(ncol(X)), S)
    Xp <- Xc[, cand, drop = FALSE]
    if (length(S)) {
      Q <- qr.Q(qr(Xc[, S, drop = FALSE]))
      Xp <- Xp - Q %*% crossprod(Q, Xp)
    }
    along <- drop(crossprod(Xp, r))
    score <- ifelse(along > 0, along^2 / pmax(colSums(Xp^2), 1e-300), -Inf)
    if (!any(is.finite(score))) break
    S <- c(S, cand[which.max(score)])
    w <- .ctCoresetNnls(Xc[, S, drop = FALSE], yc)
    r <- yc - drop(Xc[, S, drop = FALSE] %*% w)
  }
  S
}

# The coreset whose expected round cost (`.ctCoresetRoundCost()`) is least,
# refitted on every row, with that cost and its failure rate as attributes;
# NULL when none is expected to cost less than a full round, or there are too
# few held-out draws near the top to judge by. Sizes are tried along one
# forward-selection path, which grows until a size's share of the round alone
# costs more than the best found: nothing larger can be cheaper. Two sets are
# held out: the newest of the rows the sampler's rounds checked in full
# (`focus`), which are where the weight is moving to, and every fourth of the
# rest. Older checked rows are fitted like any other: held out for good, a
# coreset fitted mostly on the first rounds' tails could never reach them.
# Held out deterministically, so that the sampler's random stream is the same
# with a coreset as without.
.ctCoresetChoose <- function(X, y, rows, round, start = integer(),
  focus = integer()) {
  focus <- intersect(focus, rows)
  newest <- utils::tail(focus, round$checks)
  plain <- setdiff(rows, newest)
  held <- plain[seq_along(plain) %% 4L == 0L]
  fitrows <- setdiff(plain, held)
  rel <- numeric(length(y))
  rel[rows] <- .ctCoresetRelevance(y[rows])
  if (sum(rel[held]) < 2 || length(fitrows) < 10L) return(NULL)
  N <- ncol(X)
  fullround <- round$values + 2 * round$grads
  failrate <- function(cf, idx) {
    if (!length(idx)) return(0)
    .ctCoresetFailRate(.ctCoresetPredict(cf, X[idx, , drop = FALSE]) - y[idx], rel[idx],
      round$checks)
  }
  best <- NULL
  path <- start
  M <- max(1L, length(start))
  repeat {
    floorcost <- M / N * fullround + round$checks
    if (M > N || floorcost >= fullround || (!is.null(best) && floorcost >= best$cost)) break
    if (M > length(path)) path <- .ctCoresetPath(X, y, M, fitrows, path, rel[fitrows])
    if (M > length(path)) break
    cf <- .ctCoresetWeights(X, y, path[seq_len(M)], fitrows, rel[fitrows])
    fail <- max(failrate(cf, held), failrate(cf, newest))
    cost <- .ctCoresetRoundCost(sum(cf$w > 0), N, round, fail)
    if (is.null(best) || cost < best$cost) best <- list(M = M, cost = cost, fail = fail)
    # Steps of about a tenth of the size: fine where sizes are small, and a
    # path of hundreds of units is not walked one at a time.
    M <- M + max(1L, M %/% 10L)
  }
  if (is.null(best) || best$cost >= fullround) return(NULL)
  S <- .ctCoresetPath(X, y, best$M, rows, start, rel[rows])
  cf <- .ctCoresetWeights(X, y, S, rows, rel[rows])
  cf$cost <- best$cost
  cf$fail <- best$fail
  cf
}

# The sampler's density with a coreset behind it: `lpg` (`.ctBackendLpgFunc()`)
# with 'batch', 'gradbatch', 'cost' (full evaluations spent, a gradient as
# two), 'finalize' (the end check; see `imis_is()`) and 'coreset' (the record
# `$details$importance_sampling$coreset` reports). `lpg` itself comes back
# where a coreset cannot apply: off the Laplace route, or `control$coreset`
# FALSE.
.ctCoresetLpg <- function(fit, lpg, control = list()) {
  if (isFALSE(control$coreset) || !.ctBackendIsLaplace(fit) ||
    !is.function(attr(lpg, 'batch')) || !is.function(attr(lpg, 'gradbatch'))) return(lpg)
  tol <- .ctCoresetTol
  st <- new.env(parent = emptyenv())
  st$phase <- 'training'
  st$X <- NULL
  st$y <- numeric(0)
  st$focus <- integer(0)
  st$cost <- 0
  st$sizes <- integer(0)
  st$checks <- numeric(0)
  st$final <- NULL
  # What a round is made of: its draws (the latest round's), the gradients the
  # sampler asked for between rounds, and the checks a coreset round adds.
  st$values <- 0
  st$grads <- numeric(0)
  st$gradsnow <- 0
  st$batches <- 0L
  # The current coreset's own record: rounds run on it and checks failed.
  st$runs <- 0L
  st$fails <- 0L
  st$switches <- 0L
  st$runs_total <- 0L
  st$fails_total <- 0L
  engine <- local({
    built <- NULL
    function() {
      if (is.null(built)) built <<- list(module = .ctJuliaModule(fit$model_spec$project),
        objective = .ctJuliaObjective(fit, hybrid = TRUE))
      built
    }
  })
  terms <- function(X) {
    res <- .ctBackendJuliaValue(engine()$module$ctsem_coreset_terms(engine()$objective,
      t(X)))
    list(value = as.numeric(res$value), total = as.numeric(res$total),
      L = t(matrix(as.numeric(res$unit_loglik), ncol = nrow(X))))
  }
  learn <- function(t0, focus = FALSE) {
    ok <- is.finite(t0$total) & rowSums(!is.finite(t0$L)) == 0
    if (any(ok)) {
      if (focus) st$focus <- c(st$focus, length(st$y) + seq_len(sum(ok)))
      st$X <- rbind(st$X, t0$L[ok, , drop = FALSE])
      st$y <- c(st$y, t0$total[ok])
    }
  }
  active <- function() {
    keep <- st$cf$w > 0
    list(units = st$cf$S[keep], weights = st$cf$w[keep])
  }
  roundinfo <- function() list(values = st$values,
    grads = if (length(st$grads)) mean(st$grads) else 0.3 * st$values,
    checks = min(8, st$values))
  choose <- function(start = integer()) {
    cf <- .ctCoresetChoose(st$X, st$y, seq_along(st$y), roundinfo(), start, st$focus)
    if (is.null(cf)) return(FALSE)
    st$cf <- cf
    st$sizes <- c(st$sizes, sum(cf$w > 0))
    st$runs <- 0L
    st$fails <- 0L
    TRUE
  }
  # The current coreset's expected round cost with its own record mixed in:
  # its held-out failure rate counts as two rounds' worth of evidence.
  dearer <- function() {
    fail <- (st$fails + 2 * st$cf$fail) / (st$runs + 2)
    r <- roundinfo()
    .ctCoresetRoundCost(sum(st$cf$w > 0), ncol(st$X), r, fail) >= r$values + 2 * r$grads
  }
  coreset_value <- function(X) {
    a <- active()
    v <- as.numeric(.ctBackendJuliaValue(engine()$module$ctsem_coreset_batch(
      engine()$objective, t(X), .ctJuliaVector(as.integer(a$units)),
      .ctJuliaVector(as.numeric(a$weights))))) + st$cf$c
    st$cost <- st$cost + nrow(X) * length(a$units) / ncol(st$X)
    v
  }
  full <- function(X) {
    t0 <- terms(X)
    st$cost <- st$cost + nrow(X)
    t0
  }
  batch <- function(draws) {
    X <- as.matrix(draws)
    # A round's size, and the gradients asked for since the last one (none
    # before the first: until a round has passed, a round's gradients are
    # taken as 30% of its draws, about what 'checked''s Newton components
    # ask for).
    st$values <- nrow(X)
    if (st$batches > 0L) st$grads <- c(st$grads, st$gradsnow)
    st$batches <- st$batches + 1L
    st$gradsnow <- 0
    if (st$phase == 'training') {
      # Full rounds are what the sampler pays without a coreset, so learning
      # from them costs nothing; each one, try again for a coreset that pays.
      t0 <- full(X)
      learn(t0)
      if (!is.null(st$X) && choose(if (!is.null(st$cf)) st$cf$S else integer())) {
        st$phase <- 'coreset'
        st$switches <- st$switches + 1L
      }
      v <- t0$value
      v[!is.finite(v)] <- -1e100
      return(v)
    }
    v <- coreset_value(X)
    # A round's check: the draws the coreset ranks highest, where the weight
    # will sit, and an even spread of the rest.
    k <- min(8L, nrow(X))
    top <- order(-v)[seq_len(ceiling(k / 2))]
    rest <- setdiff(seq_len(nrow(X)), top)
    pick <- c(top, rest[unique(round(seq(1, length(rest), length.out = floor(k / 2))))])
    pick <- pick[!is.na(pick)]
    t0 <- full(X[pick, , drop = FALSE])
    good <- is.finite(t0$value) & is.finite(v[pick])
    if (any(good)) {
      err <- (v[pick] - t0$value)[good]
      st$checks <- c(st$checks, err)
      v[pick[good]] <- t0$value[good]
      learn(t0, focus = TRUE)
      rel <- .ctCoresetRelevance(t0$total[good], max(st$y))
      st$runs <- st$runs + 1L
      st$runs_total <- st$runs_total + 1L
      if (!.ctCoresetPasses(err, rel, tol$drift)) {
        # Failed. The round is evaluated in full: the sampler keeps every
        # value it is handed, and a coreset that overrated the draws it
        # ranked highest by 20 log units (a first, too small, coreset on the
        # 250-subject strong shell) misled its adaptation for thirty rounds
        # and left 1,900 draws for the end check to correct. The full round is
        # training data where the sampler now is, and the coreset is chosen
        # again from it; if none pays, full rounds until one does.
        st$fails <- st$fails + 1L
        st$fails_total <- st$fails_total + 1L
        others <- setdiff(seq_len(nrow(X)), pick)
        if (length(others)) {
          t1 <- full(X[others, , drop = FALSE])
          v[others] <- t1$value
          learn(t1)
        }
        if (!choose(st$cf$S)) st$phase <- 'training'
      } else {
        keep <- st$cf[c('cost', 'fail')]
        st$cf <- c(.ctCoresetWeights(st$X, st$y, st$cf$S, seq_along(st$y)), keep)
        # Passing, but its record may still say a full round is cheaper.
        if (dearer()) st$phase <- 'training'
      }
    }
    v[!is.finite(v)] <- -1e100
    v
  }
  gradbatch <- function(draws) {
    X <- as.matrix(draws)
    st$gradsnow <- st$gradsnow + nrow(X)
    if (st$phase != 'coreset') {
      st$cost <- st$cost + 2 * nrow(X)
      return(attr(lpg, 'gradbatch')(X))
    }
    a <- active()
    out <- matrix(as.numeric(.ctBackendJuliaValue(
      engine()$module$ctsem_coreset_batch_gradient(engine()$objective, t(X),
        .ctJuliaVector(as.integer(a$units)), .ctJuliaVector(as.numeric(a$weights))))),
      nrow = ncol(X) + 1L)
    st$cost <- st$cost + 2 * nrow(X) * length(a$units) / ncol(st$X)
    value <- out[1L, ] + st$cf$c
    value[!is.finite(value)] <- -1e100
    list(value = value, gradient = t(out[-1L, , drop = FALSE]))
  }
  # The end check. The full density where the weights sit; if the coreset was
  # off there by more than `tol$final` either side of its typical offset (a
  # constant cancels), every draw carrying the weight is corrected to the full
  # density, so the weights handed back are exact where they matter.
  finalize <- function(samples, logp, logw) {
    if (is.null(st$cf)) return(NULL)
    o <- order(-logw)
    top <- o[seq_len(min(20L, length(o)))]
    vt <- as.numeric(.ctBackendJuliaValue(engine()$module$ctsem_evaluate_batch(
      engine()$objective, t(samples[top, , drop = FALSE]))))
    st$cost <- st$cost + length(top)
    delta <- vt - logp[top]
    ok <- is.finite(delta)
    offset <- if (any(ok)) stats::median(delta[ok]) else 0
    worst <- if (any(ok)) max(abs(delta[ok] - offset)) else Inf
    st$final <- list(checked = length(top), worst = worst, corrected = 0L)
    if (worst <= tol$final) return(NULL)
    cw <- cumsum(exp(logw[o] - max(logw))) / sum(exp(logw - max(logw)))
    need <- union(o[seq_len(which(cw >= 0.999)[1L])], which(logw > max(logw) + log(1e-3)))
    need <- setdiff(need, top)
    full_all <- rep(NA_real_, length(logw))
    full_all[top] <- vt
    if (length(need)) {
      full_all[need] <- as.numeric(.ctBackendJuliaValue(engine()$module$ctsem_evaluate_batch(
        engine()$objective, t(samples[need, , drop = FALSE]))))
      st$cost <- st$cost + length(need)
    }
    d <- full_all - logp
    d[!is.finite(d)] <- offset
    st$final$corrected <- length(need) + length(top)
    logw + d
  }
  fn <- function(parm) batch(matrix(parm, 1L))
  attr(fn, 'batch') <- batch
  attr(fn, 'gradbatch') <- gradbatch
  attr(fn, 'cost') <- function() st$cost
  attr(fn, 'finalize') <- finalize
  attr(fn, 'coreset') <- function() list(
    used = !is.null(st$cf),
    units = if (!is.null(st$X)) ncol(st$X),
    size = if (!is.null(st$cf)) sum(st$cf$w > 0),
    sizes = st$sizes,
    phase = st$phase,
    switches = st$switches,
    rounds_on_coreset = st$runs_total,
    failed_checks = st$fails_total,
    expected_round_cost = if (!is.null(st$cf)) st$cf$cost,
    training_draws = length(st$y),
    check_rms = if (length(st$checks)) sqrt(mean(st$checks^2)) else NA_real_,
    check_max = if (length(st$checks)) max(abs(st$checks)) else NA_real_,
    final = st$final,
    cost = st$cost)
  fn
}
