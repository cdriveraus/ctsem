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
# by the expected cost of the sampler's work in full evaluations, unit terms
# counted: never in seconds, which other work on the machine moves, so that
# the same fit makes the same choices anywhere. It evaluates in full, as the sampler would without it, and learns from that
# until some coreset is expected to make a round cheaper; then gradients run
# on that coreset's units only, and each batch of draws does too where that is
# cheaper for its size (`ctsem_coreset_batch`, several draws at once on the
# pool a full evaluation uses): a few of its draws are checked in full, and
# the rest of the batch is evaluated in full when the check fails. A batch too
# small to pay for its check is evaluated in full, and the coreset's
# prediction at its draws, read off the unit terms, checks it for nothing.
# The size is the one whose expected round cost is least
# (`.ctCoresetChoose()`): its share of the round's values and gradients, its
# checks, and the chance a check fails times the full batch that costs. That
# chance is read from held-out errors and then from the checks the coreset
# actually fails; when they make it dearer than a full round it goes back to
# full rounds and tries again with what they teach. At the end it checks the
# draws that carry the weight and corrects them where needed, so what it hands
# back is the full density's.

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

# What a batch of `round$values` draws is expected to cost on a coreset of `M`
# of `N` units, and in full, in full evaluations. On the coreset: its share
# `M / N`, the `checks` evaluated in full, and -- with probability `fail` --
# the rest evaluated in full after a failed check. Measured on dev3, a
# coreset's draw took M/N of a full one's time to within 5% at 100, 250 and
# 1,000 units, its draws evaluated several at once.
.ctCoresetBatchCost <- function(M, N, round, fail) {
  n <- round$values
  k <- min(round$checks, n)
  c(coreset = M / N * n + k + fail * (n - k), full = n)
}

# What a sampler round -- a batch of `values` draws and `grads` gradients,
# each as two values -- is expected to cost with a coreset of `M` of `N` units:
# the batch whichever way is cheaper (`.ctCoresetBatchCost()`), and the
# gradients at its share where it has passed its last check and in full where
# it failed it, since a coreset serves gradients only where it has just
# passed one. A round in full costs `.ctCoresetFullRound()`.
.ctCoresetRoundCost <- function(M, N, round, fail) {
  share <- M / N
  2 * round$grads * (share + fail * (1 - share)) +
    min(.ctCoresetBatchCost(M, N, round, fail))
}
.ctCoresetFullRound <- function(round) round$values + 2 * round$grads

# How often a round's check would fail, from a coreset's errors `err` and
# predictions `pred` at held-out draws of relevance `rel`. The check is
# simulated as the sampler makes it: a round of `info$values` draws resampled
# from these, its `info$checks` checks the half the coreset ranks highest and
# an even spread of the rest. Ranking by the coreset's own value is what finds
# a coreset too small: it overrates the draws it puts on top, and a check
# drawn by relevance alone missed that, chose sizes that failed in the
# sampler, and cost a full round each time. On a private random stream, so
# that the sampler's own is the same with a coreset as without.
.ctCoresetFailRate <- function(err, rel, pred, info) {
  n <- length(err)
  if (!n || sum(rel) <= 0) return(0)
  had <- exists('.Random.seed', envir = globalenv(), inherits = FALSE)
  if (had) old <- get('.Random.seed', envir = globalenv(), inherits = FALSE)
  on.exit(if (had) assign('.Random.seed', old, envir = globalenv()) else
    if (exists('.Random.seed', envir = globalenv(), inherits = FALSE))
      rm('.Random.seed', envir = globalenv()))
  set.seed(20261021L)
  k <- max(1L, as.integer(min(info$checks, info$values)))
  ntop <- ceiling(k / 2)
  mean(vapply(seq_len(200L), function(b) {
    i <- sample.int(n, max(k, as.integer(info$values)), replace = TRUE)
    o <- i[order(-pred[i])]
    others <- o[-seq_len(ntop)]
    spread <- others[unique(round(seq(1, length(others), length.out = k - ntop)))]
    pick <- c(o[seq_len(ntop)], spread[!is.na(spread)])
    !.ctCoresetPasses(err[pick], rel[pick], .ctCoresetTol$drift)
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
# refitted on every row, with that cost and its failure rate as fields. When
# there is none, an empty list whose 'why' says so: too few rows to fit
# ('held'), none expected to cost less than a full round ('cost', with the
# best found), or weights all zero ('empty'). Sizes are tried along one
# forward-selection path, which grows until what a size costs with checks that
# never fail is more than the best found: nothing larger can be cheaper. Two
# sets are held out: the newest of the rows the sampler's rounds checked in full
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
  if (length(fitrows) < 10L) {
    return(structure(list(), why = 'held', relheld = sum(rel[held])))
  }
  N <- ncol(X)
  fullround <- .ctCoresetFullRound(round)
  # The failure rate a held-out set shows, with its uncertainty: the set's
  # relevance, how many draws near the top it amounts to, counts as that many
  # observations, under a uniform prior. A sampler's first rounds are tails,
  # with a draw or two near the top between them; a bar on that count kept
  # the first coreset back for two to three rounds of 100, a third of a run's
  # full evaluations, where pricing an uncertain rate lets the checks decide.
  failrate <- function(cf, idx) {
    if (!length(idx)) return(NA_real_)
    pred <- .ctCoresetPredict(cf, X[idx, , drop = FALSE])
    m <- sum(rel[idx])
    (m * .ctCoresetFailRate(pred - y[idx], rel[idx], pred, round) + 1) / (m + 2)
  }
  best <- NULL
  path <- start
  M <- max(1L, length(start))
  repeat {
    # What size M costs if its checks never fail.
    floorcost <- .ctCoresetRoundCost(M, N, round, 0)
    if (M > N || floorcost >= fullround || (!is.null(best) && floorcost >= best$cost)) break
    if (M > length(path)) path <- .ctCoresetPath(X, y, M, fitrows, path, rel[fitrows])
    if (M > length(path)) break
    cf <- .ctCoresetWeights(X, y, path[seq_len(M)], fitrows, rel[fitrows])
    # Weights that all come out zero leave a constant and no unit to
    # evaluate, which costs nothing by the formula: priced, such a coreset was
    # chosen, and the engine was then handed an empty unit list.
    if (any(cf$w > 0)) {
      fail <- max(failrate(cf, held), failrate(cf, newest), na.rm = TRUE)
      cost <- .ctCoresetRoundCost(sum(cf$w > 0), N, round, fail)
      if (is.null(best) || cost < best$cost) best <- list(M = M, cost = cost, fail = fail)
    }
    # Steps of about a tenth of the size: fine where sizes are small, and a
    # path of hundreds of units is not walked one at a time.
    M <- M + max(1L, M %/% 10L)
  }
  if (is.null(best) || best$cost >= fullround) {
    return(structure(list(), why = 'cost', best = if (!is.null(best)) best,
      fullround = fullround))
  }
  S <- .ctCoresetPath(X, y, best$M, rows, start, rel[rows])
  cf <- .ctCoresetWeights(X, y, S, rows, rel[rows])
  if (!any(cf$w > 0)) return(structure(list(), why = 'empty'))
  # The coreset is the units it was priced by, those with weight: refits on
  # later rows re-weight them, and add none. Keeping the zero-weight units of
  # the path let refits give them weight, and coresets priced at 28 units of
  # 1,000 grew to 46.
  keep <- cf$w > 0
  cf$S <- cf$S[keep]
  cf$w <- cf$w[keep]
  cf$cost <- best$cost
  cf$fail <- best$fail
  cf
}

# The sampler's density with a coreset behind it: `lpg` (`.ctBackendLpgFunc()`)
# with 'batch', 'gradbatch', 'cost' (full evaluations spent, a gradient as
# two), 'finalize' (the end check; see `imis_is()`) and 'coreset' (the record
# `$details$importance_sampling$coreset` reports). `lpg` itself comes back
# where a coreset cannot apply: off the Laplace route, `control$coreset`
# FALSE, or a density that is not the fit's own sum over units (its
# 'unitsum' mark). The coreset's terms are the fit's Laplace terms, so on any
# other density -- the quadrature continuation's -- it would sample the wrong
# one without a word.
.ctCoresetLpg <- function(fit, lpg, control = list()) {
  if (isFALSE(control$coreset) || !.ctBackendIsLaplace(fit) ||
    !isTRUE(attr(lpg, 'unitsum')) || !is.function(attr(lpg, 'batch')) ||
    !is.function(attr(lpg, 'gradbatch'))) return(lpg)
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
  # What a round is made of: a batch of draws and the gradients the sampler
  # asks for before the next, each smoothed over the batches so far. The
  # profile paths' batches are a few draws each, the rounds after them
  # hundreds, and the smoothing follows the change within a few.
  st$values <- NA_real_
  st$grads <- NA_real_
  st$gradsnow <- 0
  st$batches <- 0L
  # The current coreset's own record: rounds run on it and checks failed.
  st$runs <- 0L
  st$fails <- 0L
  st$switches <- 0L
  st$runs_total <- 0L
  st$fails_total <- 0L
  st$free_total <- 0L
  # Whether the current coreset may serve gradients unchecked, and whether
  # the last check, of whichever coreset, passed. A coreset chosen while the
  # last check passed, or before any, is trusted on its held-out errors; one
  # chosen after a failed check, where the sampler has moved past the rows it
  # was fitted on, has its gradient batches checked until one passes
  # (`gradbatch`). At 'is''s profile-path points every check failed, by 10-140
  # log units, and the paths took their gradients from it. Taking the
  # gradients of every untrusted coreset in full instead cost 'checked' a
  # round of them after each failure, 25% more evaluations at 1,000 units,
  # where the round the failure had just taught it was accurate. Choices
  # made, how many then passed, and the training rows at the last one
  # (`worth_choosing()`).
  st$trusted <- FALSE
  st$last_passed <- TRUE
  st$nchoose <- 0L
  st$proved <- 0L
  st$rows_at_choose <- 0L
  # Where the time went, for the record only: seconds and draws by kind of
  # call. Nothing is decided by them.
  st$secs <- c(full = 0, coreset = 0, gradient_full = 0, gradient_coreset = 0,
    choose = 0)
  st$draws <- c(full = 0, coreset = 0, gradient_full = 0, gradient_coreset = 0)
  # What the full draws were for: rounds while there is no coreset, batches
  # too small for one, checks of value and gradient batches, the rest of a
  # batch whose check failed, and the end check.
  st$why <- c(training = 0, small = 0, check = 0, gradient_check = 0,
    failed_rest = 0, end = 0)
  # The full values already known, by draw, so that the end check evaluates
  # only draws whose values came from a coreset. Keyed to 12 significant
  # figures: the sampler's draws can reach the end check through another
  # matrix product than the one that made them.
  st$known <- new.env(hash = TRUE, parent = emptyenv())
  st$trace <- list()
  drawkey <- function(X) apply(signif(as.matrix(X), 12), 1L, paste, collapse = ',')
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
  timed <- function(kind, expr, draws) {
    started <- proc.time()[["elapsed"]]
    out <- expr
    tally(kind, proc.time()[["elapsed"]] - started, draws)
    out
  }
  tally <- function(kind, secs, draws = 0) {
    st$secs[[kind]] <- st$secs[[kind]] + secs
    if (kind %in% names(st$draws)) st$draws[[kind]] <- st$draws[[kind]] + draws
  }
  roundinfo <- function() list(values = st$values,
    grads = if (is.na(st$grads)) 0.3 * st$values else st$grads, checks = 8)
  choose <- function(start = integer()) {
    if (!worth_choosing()) return(FALSE)
    started <- proc.time()[["elapsed"]]
    cf <- .ctCoresetChoose(st$X, st$y, seq_along(st$y), roundinfo(), start, st$focus)
    tally('choose', proc.time()[["elapsed"]] - started)
    st$nchoose <- st$nchoose + 1L
    st$rows_at_choose <- length(st$y)
    st$trace[[length(st$trace) + 1L]] <- c(rows = length(st$y),
      M = if (length(cf)) sum(cf$w > 0) else 0,
      why = match(.ctJuliaOr(attr(cf, 'why'), 'chosen'), c('chosen', 'held', 'cost', 'empty')),
      fail = if (length(cf)) cf$fail else .ctJuliaOr(attr(cf, 'best')$fail, NA),
      cost = if (length(cf)) cf$cost else .ctJuliaOr(attr(cf, 'best')$cost, NA),
      full = .ctJuliaOr(attr(cf, 'fullround'), NA),
      bestM = .ctJuliaOr(attr(cf, 'best')$M, NA),
      relheld = .ctJuliaOr(attr(cf, 'relheld'), NA))
    if (!length(cf)) return(FALSE)
    st$cf <- cf
    st$sizes <- c(st$sizes, sum(cf$w > 0))
    st$runs <- 0L
    st$fails <- 0L
    st$trusted <- st$last_passed
    TRUE
  }
  # A choice refits the forward selection on every training row, R arithmetic
  # over rows and units. It is made again only once the sampler has added at
  # least as many rows as the last coreset had units, one for each weight it
  # fits: with fewer, it refits to nearly the coreset that just failed. On the
  # profile paths of 'is', where every check failed, a choice after each took
  # a third of the run. Rows, not seconds, so the choice does not depend on
  # the machine's load.
  worth_choosing <- function() {
    if (is.null(st$cf)) return(TRUE)
    length(st$y) - st$rows_at_choose >= sum(st$cf$w > 0)
  }
  # How often the current coreset's checks fail: its held-out rate counts as
  # two checks' worth of evidence, and then its own record. Not the sampler's
  # recent checks, which may have been of other coresets elsewhere: on 'is'
  # those of the profile paths, all failed, kept every later round off a
  # coreset that would have passed.
  failnow <- function() (st$fails + 2 * st$cf$fail) / (st$runs + 2)
  # Whether the current coreset's expected round cost is now above a full
  # round's.
  dearer <- function() {
    r <- roundinfo()
    .ctCoresetRoundCost(sum(st$cf$w > 0), ncol(st$X), r, failnow()) >=
      .ctCoresetFullRound(r)
  }
  smooth <- function(old, new) if (is.na(old)) new else 0.7 * old + 0.3 * new
  coreset_value <- function(X) {
    a <- active()
    share <- length(a$units) / ncol(st$X)
    v <- timed('coreset', as.numeric(.ctBackendJuliaValue(
      engine()$module$ctsem_coreset_batch(engine()$objective, t(X),
        .ctJuliaVector(as.integer(a$units)), .ctJuliaVector(as.numeric(a$weights))))),
      nrow(X))
    st$cost <- st$cost + nrow(X) * share
    v + st$cf$c
  }
  full <- function(X, why) {
    engine()
    t0 <- timed('full', terms(X), nrow(X))
    st$cost <- st$cost + nrow(X)
    st$why[[why]] <- st$why[[why]] + nrow(X)
    ok <- is.finite(t0$value)
    if (any(ok)) list2env(stats::setNames(as.list(t0$value[ok]), drawkey(X[ok, , drop = FALSE])),
      envir = st$known)
    t0
  }
  # The full density at draws `X`: known values where a draw was evaluated in
  # full already, the rest evaluated now and counted for the end check.
  fullvalues <- function(X) {
    out <- unlist(mget(drawkey(X), envir = st$known, ifnotfound = list(NA_real_)),
      use.names = FALSE)
    todo <- which(is.na(out))
    if (length(todo)) {
      started <- proc.time()[["elapsed"]]
      out[todo] <- as.numeric(.ctBackendJuliaValue(engine()$module$ctsem_evaluate_batch(
        engine()$objective, t(X[todo, , drop = FALSE]))))
      st$cost <- st$cost + length(todo)
      st$why[['end']] <- st$why[['end']] + length(todo)
      tally('full', proc.time()[["elapsed"]] - started, length(todo))
    }
    out
  }
  part <- function(t0, i) list(value = t0$value[i], total = t0$total[i],
    L = t0$L[i, , drop = FALSE])
  # A check: the draws the coreset ranks highest by `v`, where the weight
  # will sit, and an even spread of the rest.
  checkpick <- function(v) {
    k <- min(8L, length(v))
    top <- order(-v)[seq_len(ceiling(k / 2))]
    rest <- setdiff(seq_along(v), top)
    pick <- c(top, rest[unique(round(seq(1, length(rest), length.out = floor(k / 2))))])
    pick[!is.na(pick)]
  }
  # A failed check: the coreset is chosen again from what the sampler has
  # taught it, and if none pays, full rounds until one does.
  failed <- function() {
    st$fails <- st$fails + 1L
    st$fails_total <- st$fails_total + 1L
    st$trusted <- FALSE
    st$last_passed <- FALSE
    if (!choose(st$cf$S)) st$phase <- 'training'
  }
  passed <- function() {
    if (!st$trusted) st$proved <- st$proved + 1L
    st$trusted <- TRUE
    st$last_passed <- TRUE
    keep <- st$cf[c('cost', 'fail')]
    refit <- .ctCoresetWeights(st$X, st$y, st$cf$S, seq_along(st$y))
    # Refitted on the new rows, its weights can all come out zero: chosen
    # again, as after a failure.
    if (!any(refit$w > 0)) {
      if (!choose(st$cf$S)) st$phase <- 'training'
      return(invisible())
    }
    st$cf <- c(refit, keep)
    # Passing, but its record may still say a full round is cheaper.
    if (dearer()) st$phase <- 'training'
  }
  batch <- function(draws) {
    X <- as.matrix(draws)
    n <- nrow(X)
    # A batch's size, and the gradients asked for since the last one (none
    # before the first: until then a round's gradients are taken as 30% of
    # its draws, about what 'checked''s Newton components ask for).
    st$values <- smooth(st$values, n)
    if (st$batches > 0L) st$grads <- smooth(st$grads, st$gradsnow)
    st$batches <- st$batches + 1L
    st$gradsnow <- 0
    v <- if (st$phase == 'coreset') withcoreset(X) else inparts(X, 'training')
    v[!is.finite(v)] <- -1e100
    v
  }
  # A batch while there is a coreset: on it where that pays for the batch's
  # size, else in full.
  withcoreset <- function(X) {
    n <- nrow(X)
    this <- roundinfo()
    this$values <- n
    bc <- .ctCoresetBatchCost(sum(st$cf$w > 0), ncol(st$X), this, failnow())
    if (bc[['coreset']] < bc[['full']]) return(oncoreset(X))
    t0 <- full(X, 'small')
    # In full because the coreset would not pay for its check on a batch this
    # small, so the check is free: the coreset's prediction from the unit
    # terms, at the draws a check would read.
    pred <- .ctCoresetPredict(st$cf, t0$L)
    pick <- checkpick(ifelse(is.finite(pred), pred, -Inf))
    good <- is.finite(t0$value[pick]) & is.finite(pred[pick])
    learn(part(t0, pick), focus = TRUE)
    learn(part(t0, setdiff(seq_len(n), pick)))
    if (any(good)) {
      err <- (pred[pick] - t0$value[pick])[good]
      st$checks <- c(st$checks, err)
      st$runs <- st$runs + 1L
      st$free_total <- st$free_total + 1L
      rel <- .ctCoresetRelevance(t0$total[pick][good], max(st$y))
      if (!.ctCoresetPasses(err, rel, tol$drift)) failed() else passed()
    }
    t0$value
  }
  # A batch with no coreset to run on -- the first rounds, or the rest of one
  # whose check failed -- is evaluated in full a half at a time, while the
  # other half could pay for a check on a coreset (on one of the last
  # coreset's size that never fails). After each half a coreset is chosen
  # from what the sampler has taught so far, and the rest runs on it if one
  # pays. Full rounds are what the sampler pays without a coreset, so learning
  # from them costs nothing; evaluated whole, they and failed batches were a
  # third to a half of a run's full evaluations.
  inparts <- function(X, why) {
    n <- nrow(X)
    v <- rep(NA_real_, n)
    done <- 0L
    while (done < n) {
      left <- n - done
      half <- left %/% 2L
      r <- roundinfo()
      r$values <- half
      bc <- if (is.null(st$X)) .ctCoresetBatchCost(0, 1, r, 0) else
        .ctCoresetBatchCost(if (is.null(st$cf)) 0 else sum(st$cf$w > 0), ncol(st$X), r, 0)
      take <- if (half > 0L && bc[['coreset']] < bc[['full']]) left - half else left
      idx <- done + seq_len(take)
      t0 <- full(X[idx, , drop = FALSE], why)
      v[idx] <- t0$value
      learn(t0)
      done <- done + take
      if (!is.null(st$X) && choose(if (!is.null(st$cf)) st$cf$S else integer())) {
        st$phase <- 'coreset'
        st$switches <- st$switches + 1L
        if (done < n) {
          rest <- seq.int(done + 1L, n)
          v[rest] <- withcoreset(X[rest, , drop = FALSE])
        }
        return(v)
      }
    }
    v
  }
  oncoreset <- function(X) {
    v <- coreset_value(X)
    pick <- checkpick(v)
    t0 <- full(X[pick, , drop = FALSE], 'check')
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
        # Failed. None of the round's values comes from this coreset: the
        # sampler keeps every value it is handed, and a coreset that overrated
        # the draws it ranked highest by 20 log units (a first, too small,
        # coreset on the 250-subject strong shell) misled its adaptation for
        # thirty rounds and left 1,900 draws for the end check to correct. The
        # rest is evaluated as a batch with no coreset is (`inparts`): its
        # full halves are training data where the sampler now is, and a
        # coreset chosen from them takes what is left, checked again.
        st$fails <- st$fails + 1L
        st$fails_total <- st$fails_total + 1L
        st$trusted <- FALSE
        st$last_passed <- FALSE
        st$phase <- 'training'
        others <- setdiff(seq_len(nrow(X)), pick)
        if (length(others)) {
          v[others] <- inparts(X[others, , drop = FALSE], 'failed_rest')
        } else if (choose(st$cf$S)) st$phase <- 'coreset'
      } else passed()
    }
    v[!is.finite(v)] <- -1e100
    v
  }
  gradfull <- function(X) {
    started <- proc.time()[["elapsed"]]
    st$cost <- st$cost + 2 * nrow(X)
    out <- attr(lpg, 'gradbatch')(X)
    tally('gradient_full', proc.time()[["elapsed"]] - started, nrow(X))
    out
  }
  gradbatch <- function(draws) {
    X <- as.matrix(draws)
    n <- nrow(X)
    st$gradsnow <- st$gradsnow + n
    if (st$phase != 'coreset') return(gradfull(X))
    a <- active()
    share <- length(a$units) / ncol(st$X)
    # A coreset not yet trusted here has its gradient batch checked as a
    # batch of values is: its values at a check's draws evaluated in full, and
    # the batch's gradients taken in full if they fail. Where that is not
    # expected to cost less than the batch in full -- a few gradients, or a
    # coreset that keeps failing -- they are taken in full straight away.
    k <- min(8L, n)
    if (!st$trusted && 2 * share * n + k + failnow() * 2 * n >= 2 * n) return(gradfull(X))
    started <- proc.time()[["elapsed"]]
    out <- matrix(as.numeric(.ctBackendJuliaValue(
      engine()$module$ctsem_coreset_batch_gradient(engine()$objective, t(X),
        .ctJuliaVector(as.integer(a$units)), .ctJuliaVector(as.numeric(a$weights))))),
      nrow = ncol(X) + 1L)
    tally('gradient_coreset', proc.time()[["elapsed"]] - started, n)
    st$cost <- st$cost + 2 * n * share
    value <- out[1L, ] + st$cf$c
    if (!st$trusted) {
      pick <- checkpick(ifelse(is.finite(value), value, -Inf))
      t0 <- full(X[pick, , drop = FALSE], 'gradient_check')
      good <- is.finite(t0$value) & is.finite(value[pick])
      if (any(good)) {
        err <- (value[pick] - t0$value)[good]
        st$checks <- c(st$checks, err)
        learn(t0, focus = TRUE)
        st$runs <- st$runs + 1L
        rel <- .ctCoresetRelevance(t0$total[good], max(st$y))
        if (!.ctCoresetPasses(err, rel, tol$drift)) {
          failed()
          return(gradfull(X))
        }
        value[pick[good]] <- t0$value[good]
        passed()
      }
    }
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
    vt <- fullvalues(samples[top, , drop = FALSE])
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
    if (length(need)) full_all[need] <- fullvalues(samples[need, , drop = FALSE])
    d <- full_all - logp
    d[!is.finite(d)] <- offset
    st$final$corrected <- length(need) + length(top)
    logw + d
  }
  fn <- function(parm) batch(matrix(parm, 1L))
  attr(fn, 'batch') <- batch
  # The full objective's engine, for a stage that runs there on the full
  # density (`.ctDaisEngine()`); the coreset serves only the batches here.
  attr(fn, 'engine') <- attr(lpg, 'engine')
  attr(fn, 'gradbatch') <- gradbatch
  # The full density's gradient, for a check of the answer
  # (`.ctSteinCheck()`) that a coreset's error must not reach.
  attr(fn, 'fullgradbatch') <- gradfull
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
    free_checks = st$free_total,
    chosen = st$nchoose,
    chosen_passed = st$proved,
    seconds = st$secs,
    draws = st$draws,
    full_for = st$why,
    choices = if (length(st$trace)) do.call(rbind, st$trace),
    training_draws = length(st$y),
    check_rms = if (length(st$checks)) sqrt(mean(st$checks^2)) else NA_real_,
    check_max = if (length(st$checks)) max(abs(st$checks)) else NA_real_,
    final = st$final,
    cost = st$cost)
  fn
}
