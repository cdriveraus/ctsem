# Doubly adaptive importance sampling: the normal moved to the posterior's
# moments.
#
# `uncertainty = 'checked'` tests the Hessian's normal with gradients at draws
# from it (`.ctOptimNormalCheck()`), takes one step of the posterior's mean
# from the pairs' intercept and one of its widths from their odd and even
# parts, and hands the directions it names to importance sampling with the
# rest held at that normal. Where the posterior's mass is far from the mode --
# a 75-parameter time series whose bulk sits 95 normal sds from it, 41 of its
# directions three times wider than the Hessian says -- one step is not
# enough, importance sampling in 67 named directions kept 2 draws in 2,000,
# and what was reported was the Hessian's normal with a warning.
#
# This iterates that step with damping until the normal sits on the
# posterior's moments: van den Boom, Cremaschi and Thiery, "Doubly adaptive
# importance sampling", JCGS 35 (2026), arXiv:2404.18556. In the normal's
# whitened coordinates `z` (`theta = centre + A z`, `A A'` the covariance),
# with `q = N(mu, Gamma)` the current normal there, each iteration draws `S`
# points from `q` in antithetic pairs, takes their log densities and
# gradients, and forms `Phi = log p - log q`. The damping `gamma` is the
# largest value in (0, 1] whose weights `exp(gamma Phi)` keep an effective
# sample of half the draws (bisection on the same draws). Stein's lemma for
# the Gaussian `q` then gives the damped target `q exp(gamma Phi)` its moments
# from the weighted gradients,
#
#     mu'    = mu    + gamma E_w[Gamma grad Phi]
#     Gamma' = Gamma + gamma Cov_w(Gamma grad Phi, z),
#
# with `Gamma grad Phi(z) = Gamma grad log p(z) + (z - mu)`, and `q` moves to
# them. At `gamma = 1` that is the posterior's mean and covariance by
# importance sampling; at `gamma -> 0` it is natural-gradient Gaussian
# variational inference (Khan and Rue 2023's Bayesian learning rule; the fixed
# point `E_q[grad log p] = 0`, `Gamma^-1 = E_q[-Hessian]` of Opper and
# Archambeau 2009). Below one, which is where it stays in many dimensions --
# a Gaussian cannot weight a skewed 75-dimensional posterior to half its draws
# -- the fixed point is the bridge `q^(1 - gamma) p^gamma`'s moments: the mean
# the posterior's, the sds a little under it. Measured
# (review/SAMPLING-smc-2026-10-23.md, scripts/gsm-uncertainty-2026-10/
# hd_dais.jl): on a rotated log-gamma skew in 10 of 75 dimensions the Laplace
# normal's worst mean error of 0.67 sds became 0.08 and its sds 0.90-1.01 of
# the truth at 8.6k gradients; on the time series every parameter's mean came
# within 0.4 reference sds (the check's normal: 2.5) and the Stein check
# passed, at 38-77k gradients, a seventh of a sampler's, with the sds 0.76 of
# the posterior's (0.41 in the worst direction) at any damping. The sampled
# directions' importance sampling then starts from a normal near the mass,
# which halved its cost on the 266-parameter fit of
# review/UNCERTAINTY-mean-shift-2026-10-22.md; where it still fails, this
# normal is reported instead of the Hessian's.
#
# Three things the paper's 10^5 draws an iteration never met at 64-256:
#
#   * The iterate is noisy where it is not wrong: its per-coordinate sds come
#     out right while its eigenstructure is noise (an exact 100-dimensional
#     normal's iterate gave importance weights with an effective sample of 37
#     in 128). So once settled, `average` more iterations are run and their
#     iterates averaged (Polyak), and the average is the answer; the last
#     `average` are averaged where the budget ends first.
#   * When to stop. Four rules failed, each on the time series: `gamma`
#     reaching 1 (never, in 75 dimensions); the damped step being small (it is
#     whenever `gamma` is); the undamped step against a crude noise floor
#     (within it at every iteration while the cloud moved 0.3 sds a step in
#     one direction for 90 iterations); the ELBO's drift over five iterations
#     (far from the posterior one ELBO estimate is noisy enough to hide 0.4
#     nats a step). What holds is three tests together over a window of
#     `window` iterations, `settle` times running: the natural gradient
#     `E_w[Gamma grad Phi]` within twice its own weighted standard error; the
#     ELBO's last window within two standard errors of the one before; and
#     the mean's displacement over the window within twice the random walk
#     its damped steps' noise would make. A per-iteration test cannot see
#     steps each within its noise that all point one way.
#   * A covariance step that is not positive definite halves `gamma` on the
#     same draws until it is.
#
# Cost: `iterations x S` gradients, one bridge call an iteration through the
# density's `'gradbatch'` (`.ctBackendLpgFunc()`); 60-150 iterations on the
# targets measured. On the Laplace route the density is `.ctCoresetLpg()`'s,
# whose `'gradbatch'` serves gradients on a coreset of the units where it
# has one it trusts. A coreset learns only from value batches, so where the
# density carries one (`'coreset'`), each iteration's draws are evaluated as
# a value batch too: in full while the coreset trains, a third more than the
# gradients alone, and on the coreset with its checks once chosen, which is
# what keeps the coreset honest as the normal moves -- it was fitted at the
# check's draws, and on a shell this moved the normal to where no training
# draw was. The values the batch returns are the ones weighted.
#
# A direction the data do not identify has no fixed point here: the
# gradient along it is zero, so each step multiplies its variance by
# `1 + gamma`, and on a fixture with one such direction the variance reached
# 5e8 times the curvature's in 100 iterations and the Laplace correction's
# proposal, which starts from the fit's covariance, was useless. So the
# directions the flat screen confirmed (`hold`) are held at the normal's
# width -- the moves and the covariance are projected off them -- and a
# covariance that still grows past `1e4` times the start's in any direction
# ends the run as diverged, `used = FALSE`, and the check's normal stands.
#
# Nothing here runs where the check named no direction: the Hessian's normal
# is the answer there, and this costs nothing. `control$dais = FALSE` skips
# it, which is also the oracle for measuring it.

# Draws an iteration: the dimension, at least 64 and at most 256, even.
.ctDaisDraws <- function(k) as.integer(2 * ceiling(min(256, max(64, k)) / 2))

# The iteration described above, from the normal `N(centre, cov)` on the
# density `lpg` (its `'gradbatch'` attribute; NULL without one). Returns the
# settled normal in the raw coordinates (`mean`, `cov`), where the check's
# step and recheck then take over, with its record: iterations, how many were
# averaged, whether it settled or ran out, the last damping and effective
# sample, the draws an iteration, the gradients and values spent, the cost in
# full evaluations (the density's own count where it keeps one, a gradient as
# two otherwise), and the trace per iteration.
#' @keywords internal
.ctDaisRun <- function(lpg, centre, cov, control = list(), verbose = 0,
  hold = NULL, rtol = .ctFlatDirectionRtol(), maxgrowth = 1e4) {
  gradbatch <- attr(lpg, 'gradbatch')
  if (!is.function(gradbatch)) return(NULL)
  teach <- if (is.function(attr(lpg, 'coreset'))) attr(lpg, 'batch') else NULL
  costfun <- attr(lpg, 'cost')
  cost0 <- if (is.function(costfun)) costfun() else 0
  centre <- as.numeric(centre)
  cov <- as.matrix(cov)
  e <- eigen((cov + t(cov)) / 2, symmetric = TRUE)
  keep <- e$values > rtol * max(e$values)
  A <- sweep(e$vectors[, keep, drop = FALSE], 2, sqrt(e$values[keep]), '*')
  k <- ncol(A)
  if (!k) return(NULL)
  # The held directions in the whitened coordinates, orthonormal: `A^+ h`
  # for each, where `A^+ = D^-1/2 V'` on the kept eigenvectors.
  U <- NULL
  if (length(hold) && NCOL(hold)) {
    Uz <- sweep(crossprod(e$vectors[, keep, drop = FALSE], as.matrix(hold)), 1,
      sqrt(e$values[keep]), '/')
    q <- qr(Uz)
    if (q$rank > 0) U <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  }
  offhold <- function(M) if (is.null(U)) M else M - (M %*% U) %*% t(U)
  S <- as.integer(.ctJuliaOr(control$daisDraws, .ctDaisDraws(k)))
  S <- max(4L, 2L * ((S + 1L) %/% 2L))
  maxiter <- max(1L, as.integer(.ctJuliaOr(control$daisMaxIter, 150L)))
  essfrac <- 0.5
  settle <- 3L
  window <- 20L
  average <- 20L
  ess_logw <- function(lw) {
    m <- max(lw)
    if (!is.finite(m)) return(0)
    w <- exp(lw - m)
    sum(w)^2 / sum(w^2)
  }
  mu <- numeric(k)
  Gamma <- diag(k)
  elbos <- numeric()
  mus <- list()
  stepnoise <- numeric()
  recent <- list()
  averaging <- 0L
  settled <- 0L
  musum <- numeric(k)
  Gsum <- matrix(0, k, k)
  gradients <- 0L
  values <- 0L
  trace <- list()
  gamma <- NA_real_
  essw <- NA_real_
  unusable <- 0L
  for (iter in seq_len(maxiter)) {
    Lc <- t(chol(Gamma))
    half <- S %/% 2L
    E <- matrix(stats::rnorm(half * k), half, k)
    E <- rbind(E, -E)
    Z <- sweep(E %*% t(Lc), 2, mu, '+')
    lq <- -0.5 * rowSums(E^2) - sum(log(diag(Lc))) - 0.5 * k * log(2 * pi)
    X <- sweep(Z %*% t(A), 2, centre, '+')
    out <- gradbatch(X)
    gradients <- gradients + S
    G <- as.matrix(out$gradient) %*% A
    value <- as.numeric(out$value)
    if (!is.null(teach)) {
      value <- as.numeric(teach(X))
      values <- values + S
    }
    ok <- is.finite(value) & value > -1e99 & rowSums(!is.finite(G)) == 0
    unusable <- unusable + sum(!ok)
    # A draw the density cannot use carries no weight, and no gradient
    # either: zero weight times an infinite gradient is NaN, which reached the
    # weighted means and stopped a fit on the stopping rule's NA
    # (test-poprank.R's full-rank augmented fit, on Linux).
    G[!ok, ] <- 0
    if (sum(ok) < 2L) {
      if (verbose > 0) message(sprintf(
        'DAIS: %d of %d draws usable at iteration %d; the check\'s normal stands.',
        sum(ok), S, iter))
      return(list(used = FALSE, unusable = unusable, iterations = iter,
        draws = S, gradients = gradients, values = values,
        cost = if (is.function(costfun)) costfun() - cost0 else 2 * gradients + values,
        trace = if (length(trace)) do.call(rbind, trace)))
    }
    Phi <- ifelse(ok, value - lq, -Inf)
    # Off the held directions, where the normal stands as it is.
    GP <- offhold(G %*% Gamma + sweep(Z, 2, mu))
    gamma <- 1
    if (ess_logw(Phi) < essfrac * S) {
      lo <- 0
      hi <- 1
      for (b in seq_len(60)) {
        mid <- (lo + hi) / 2
        if (ess_logw(mid * Phi) >= essfrac * S) lo <- mid else hi <- mid
      }
      gamma <- max(lo, 1e-12)
    }
    for (attempt in seq_len(40)) {
      lw <- gamma * Phi
      w <- exp(lw - max(lw))
      w[!ok] <- 0
      w <- w / sum(w)
      gbar <- colSums(w * GP)
      zbar <- colSums(w * Z)
      munew <- mu + gamma * gbar
      Cv <- crossprod(w * sweep(GP, 2, gbar), offhold(sweep(Z, 2, zbar)))
      Gnew <- Gamma + gamma * (Cv + t(Cv)) / 2
      Gnew <- (Gnew + t(Gnew)) / 2
      if (!inherits(try(chol(Gnew), silent = TRUE), 'try-error')) break
      gamma <- gamma / 2
      if (attempt == 40L) {
        Gnew <- Gamma
        munew <- mu
      }
    }
    essw <- ess_logw(gamma * Phi)
    shift <- sqrt(sum(forwardsolve(Lc, munew - mu)^2))
    natural <- sqrt(sum(forwardsolve(Lc, gbar)^2))
    W <- forwardsolve(Lc, t(sweep(GP, 2, gbar)))
    se2 <- sum(w^2 * colSums(W^2))
    noise <- 2 * sqrt(se2)
    lpf <- value[ok]
    elbo <- mean(lpf) + sum(log(diag(Lc)))
    elbos <- c(elbos, elbo)
    elbo_se <- stats::sd(lpf) / sqrt(length(lpf))
    m <- window
    n_e <- length(elbos)
    elbo_settled <- n_e >= 2L * m && (mean(elbos[(n_e - m + 1L):n_e]) -
      mean(elbos[(n_e - 2L * m + 1L):(n_e - m)])) <= 2 * elbo_se * sqrt(2 / m)
    stepnoise <- c(stepnoise, gamma * sqrt(se2))
    travel_settled <- length(mus) >= m &&
      sqrt(sum(forwardsolve(Lc, munew - mus[[length(mus) - m + 1L]])^2)) <=
        2 * sqrt(sum(utils::tail(stepnoise, m)^2))
    mu <- munew
    Gamma <- Gnew
    mus[[length(mus) + 1L]] <- mu
    recent[[length(recent) + 1L]] <- list(mu = mu, Gamma = Gamma)
    if (length(recent) > average) recent <- recent[-1L]
    trace[[iter]] <- c(gamma = gamma, ess = essw, natural = natural, noise = noise,
      elbo = elbo, shift = shift, averaging = as.numeric(averaging > 0L),
      usable = sum(ok))
    # A variance still growing past `maxgrowth` times the start's is a
    # direction with no fixed point: the run ends and the caller keeps the
    # check's normal.
    growth <- max(eigen(Gamma, symmetric = TRUE, only.values = TRUE)$values)
    if (growth > maxgrowth) {
      if (verbose > 0) message(sprintf(
        'DAIS: the covariance grew %.0f times the normal\'s along a direction after %d iterations; the data do not identify it, and the check\'s normal stands.',
        growth, iter))
      return(list(used = FALSE, diverged = TRUE, iterations = iter, growth = growth,
        unusable = unusable, draws = S, gradients = gradients, values = values,
        cost = if (is.function(costfun)) costfun() - cost0 else 2 * gradients + values,
        trace = do.call(rbind, trace)))
    }
    if (verbose > 1) message(sprintf(
      'DAIS %d: gamma %.3g, ESS %.0f of %d, natural step %.2f (noise %.2f), ELBO %.2f%s',
      iter, gamma, essw, S, natural, noise, elbo, if (averaging > 0L) ', averaging' else ''))
    # `averaging` counts the iterates summed, and is the divisor: counted as
    # one more than the terms, it inflated every averaged shift and covariance
    # by 21/20, which a held direction's variance of 1.05 showed.
    if (averaging == 0L) {
      # Each test read as unmet where it cannot be read.
      settled <- if (isTRUE(natural <= noise) && isTRUE(elbo_settled) &&
        isTRUE(travel_settled)) settled + 1L else 0L
    }
    if (averaging > 0L || settled >= settle) {
      averaging <- averaging + 1L
      musum <- musum + mu
      Gsum <- Gsum + Gamma
      if (averaging >= average) break
    }
  }
  if (averaging > 0L) {
    mu <- musum / averaging
    Gamma <- Gsum / averaging
  } else if (length(recent)) {
    mu <- Reduce(`+`, lapply(recent, `[[`, 'mu')) / length(recent)
    Gamma <- Reduce(`+`, lapply(recent, `[[`, 'Gamma')) / length(recent)
  }
  Gamma <- (Gamma + t(Gamma)) / 2
  covraw <- A %*% Gamma %*% t(A)
  cost <- if (is.function(costfun)) costfun() - cost0 else 2 * gradients + values
  if (verbose > 0) message(sprintf(
    'DAIS: %d iterations (%d averaged), %d gradients, last damping %.2g; the normal %s.',
    length(trace), averaging, gradients, gamma,
    if (averaging > 0L) 'settled' else 'ran out of iterations'))
  list(used = TRUE, mean = centre + as.numeric(A %*% mu), cov = (covraw + t(covraw)) / 2,
    iterations = length(trace), averaged = averaging, settled = averaging > 0L,
    gamma = gamma, ess = essw, unusable = unusable, draws = S, gradients = gradients,
    values = values, cost = cost, trace = do.call(rbind, trace))
}
