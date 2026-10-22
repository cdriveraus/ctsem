# Hamiltonian Monte Carlo with many chains in lockstep, for the directions
# `'checked'` samples when importance sampling cannot: too many of them, or a
# proposal its weights reject. Every leapfrog step is one call of a density's
# `'gradbatch'` for all chains -- one bridge call on julia -- which is how a
# ctsem density is cheap to ask; a chain at a time would pay the bridge per
# step.
#
# Why. Importance sampling's weights are a product over the sampled
# directions, so their spread grows with the number of them however well
# each one's proposal fits, and past about fifty it collapses. On 100
# independent log-gamma coordinates, 51 of them named, it stopped after 100
# rounds at an effective sample of 16, worst mean 0.45 sds off; on a single
# 400-row six-process time series (75 parameters, 55 named) at 2, Pareto k
# 1.5. This reached 214 and 268, means within 0.15 sds
# (scripts/gsm-uncertainty-2026-10/proto_hmc.R, proto_hmc_fit.R). Below that
# importance sampling is the cheaper: at 32 named directions 5,210
# evaluations against 13,026, and every evaluation here is a gradient, which
# costs 2-3.5 values on julia (time_valgrad.R).
#
# One step size and one trajectory length for every chain at a time, each
# jittered (the chains of Hoffman, Radul & Sountsov 2021 run so). Warm-up
# adapts the step size by dual averaging to an acceptance of `delta`
# (Hoffman & Gelman 2014) and, from its second half, a diagonal metric from
# the chains' variances and the trajectory length from the largest sd, pi /
# 2 of it. Sampling then runs in blocks of `block` until every coordinate's
# bulk and tail ESS (rstan's, over the chains) reaches `target_ess` and its
# R-hat is under `rhat`, or `max_iter` iterations. The coordinates are where
# `'checked'`'s normal is roughly N(0, I), so warm-up starts close.
#
# A curved posterior -- a ring, a funnel -- defeats any one metric: on ten
# rings or a ten-dimensional funnel the acceptance went to zero. A block
# accepting under 5% ends the run, which says so (`converged`), and the
# caller keeps the normal with a warning. `max_iter` iterations after
# warm-up at most: 6,400 draws at 16 chains.
#
# Returns the draws (chains' iterations stacked), the gradients spent, ESS,
# R-hat, acceptance and the settings it adapted.
.ctHmcLockstep <- function(gradbatch, k, chains = 16L, warmup = 60L, target_ess = 200,
  rhat = 1.05, max_iter = 400L, delta = 0.8, block = 10L, max_steps = 64L) {
  M <- chains
  z <- matrix(stats::rnorm(M * k), M, k)
  ev <- gradbatch(z)
  lv <- as.numeric(ev$value)
  g <- as.matrix(ev$gradient)
  ngrad <- M
  bad <- !is.finite(lv) | lv < -1e99 | rowSums(!is.finite(g)) > 0
  if (any(bad)) {
    z[bad, ] <- 0
    ev0 <- gradbatch(z[bad, , drop = FALSE])
    lv[bad] <- as.numeric(ev0$value)
    g[bad, ] <- as.matrix(ev0$gradient)
    ngrad <- ngrad + sum(bad)
  }
  mass <- rep(1, k)
  eps <- 0.5 / k^0.25
  len <- 1.5
  mu <- log(10 * eps)
  hbar <- 0
  lbar <- 0
  draws <- list()
  acc_all <- numeric(0)
  diag <- NULL
  it <- 0L
  steps <- function(e) max(1L, min(max_steps, as.integer(ceiling(len / e))))
  repeat {
    it <- it + 1L
    warm <- it <= warmup
    e_it <- eps * stats::runif(1, 0.8, 1.2)
    L <- max(1L, min(max_steps, as.integer(ceiling(stats::runif(1, 0.5, 1.5) * len / e_it))))
    p0 <- matrix(stats::rnorm(M * k), M, k) * rep(sqrt(mass), each = M)
    p <- p0 + e_it / 2 * g
    zn <- z
    for (l in seq_len(L)) {
      zn <- zn + e_it * sweep(p, 2, mass, '/')
      ev <- gradbatch(zn)
      ngrad <- ngrad + M
      gn <- as.matrix(ev$gradient)
      lvn <- as.numeric(ev$value)
      broken <- !is.finite(lvn) | lvn < -1e99 | rowSums(!is.finite(gn)) > 0
      gn[broken, ] <- 0
      if (l < L) p <- p + e_it * gn
    }
    p <- p + e_it / 2 * gn
    loga <- (lvn - rowSums(sweep(p^2, 2, mass, '/')) / 2) -
      (lv - rowSums(sweep(p0^2, 2, mass, '/')) / 2)
    loga[broken | !is.finite(loga)] <- -Inf
    acc <- exp(pmin(0, loga))
    take <- log(stats::runif(M)) < loga
    z[take, ] <- zn[take, ]
    lv[take] <- lvn[take]
    g[take, ] <- gn[take, ]
    if (warm) {
      t0 <- 10
      hbar <- (1 - 1 / (it + t0)) * hbar + (delta - mean(acc)) / (it + t0)
      leps <- mu - sqrt(it) / 0.05 * hbar
      eps <- exp(leps)
      lbar <- it^-0.75 * leps + (1 - it^-0.75) * lbar
      if (it > warmup / 2) draws[[length(draws) + 1L]] <- z
      if (it == warmup) {
        W <- do.call(rbind, draws)
        # Shrunk towards the start's unit variances, as Stan's windows are.
        v <- (nrow(W) * apply(W, 2, stats::var) + 5) / (nrow(W) + 5)
        mass <- 1 / v
        len <- pi / 2 * sqrt(max(v))
        eps <- exp(lbar)
        draws <- list()
      }
    } else {
      draws[[length(draws) + 1L]] <- z
      acc_all <- c(acc_all, mean(acc))
      if (length(draws) %% block == 0L || it >= warmup + max_iter) {
        arr <- simplify2array(draws)
        diag <- vapply(seq_len(k), function(j) {
          s <- t(arr[, j, ])
          c(rstan::ess_bulk(s), rstan::ess_tail(s), rstan::Rhat(s))
        }, numeric(3))
        diag[!is.finite(diag)] <- c(0, 0, Inf)[row(diag)[!is.finite(diag)]]
        # A block that hardly moves -- a curved posterior no one metric fits
        # -- will not move later either.
        if ((min(diag[1:2, ]) >= target_ess && max(diag[3, ]) <= rhat) ||
          it >= warmup + max_iter || mean(utils::tail(acc_all, block)) < 0.05) break
      }
    }
  }
  list(draws = do.call(rbind, draws), gradients = ngrad, ess = min(diag[1:2, ]),
    ess_bulk = diag[1, ], ess_tail = diag[2, ], rhat = max(diag[3, ]),
    acceptance = mean(acc_all), stepsize = eps, steps = steps(eps), mass = mass,
    iterations = length(draws), chains = M, converged = max(diag[3, ]) <= rhat &&
      min(diag[1:2, ]) >= target_ess)
}

# `'checked'`'s sampled directions by `.ctHmcLockstep()`: the same density
# `.ctImisRun()` weighs -- the directions whitened by their measured widths,
# the rest held at `centre` and tilted by its volume (`.ctImisVolumeTilt()`)
# -- with the normal's draws added in the rest. Returns what `.ctImisRun()`
# does, as far as its readers go: `theta` (`finishsamples` draws taken evenly
# through the chains), `mean`, `covariance` (the draws' plus the rest's),
# `ess`, `cost` (a gradient as two, the tilt's included), and the sampler's
# own record as `hmc`.
.ctHmcRun <- function(lpg, centre, cov, directions, widths, finishsamples = 1000L,
  target_ess = .ctEssTarget, rtol = .ctFlatDirectionRtol(), ...) {
  centre <- as.numeric(centre)
  cov <- as.matrix(cov)
  gradbase <- attr(lpg, 'gradbatch')
  if (!is.function(gradbase)) stop('HMC needs a gradient.', call. = FALSE)
  subspace <- list(V = directions, d = widths, k = ncol(directions))
  density <- .ctImisWhitenDensity(lpg, centre, subspace)
  volume <- .ctImisVolumeTilt(gradbase, centre, cov, subspace, rtol)
  if (!is.null(volume)) density <- .ctImisTiltDensity(density, volume$b)
  h <- .ctHmcLockstep(attr(density, 'gradbatch'), subspace$k, target_ess = target_ess, ...)
  Z <- h$draws
  pick <- unique(round(seq(1, nrow(Z), length.out = min(finishsamples, nrow(Z)))))
  rest <- cov - tcrossprod(directions)
  draws <- .ctImisUnwhitenMatrix(Z, centre, subspace)
  theta <- draws[pick, , drop = FALSE] + ctOptimNormalDraws(rep(0, length(centre)), rest,
    length(pick))
  gradients <- h$gradients + if (is.null(volume)) 0 else volume$gradients
  list(theta = theta, mean = colMeans(draws), covariance = stats::cov(draws) + rest,
    ess = h$ess, cost = 2 * gradients, evaluations = 0, gradients = gradients,
    tilt = if (!is.null(volume)) volume[c('size', 'noise', 'gradients')],
    hmc = h[c('ess_bulk', 'ess_tail', 'rhat', 'acceptance', 'stepsize', 'steps',
      'iterations', 'chains', 'converged')])
}
