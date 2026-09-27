# The bench's problems: a simulator and a model builder per `model` key.
#
# Every simulator is its own code with its own seed and none calls ctGenerate,
# whose draw stream moves under unrelated commits. Data are simulated once and
# stored (see `bench_data()`), so every build a cell is run against sees the
# same rows even if R's RNG or a helper used here were to change.
#
# Sources, so a number here can be traced to the study that set it:
#   stochopt regimes      dev/stochopt/models.R (copied, same draws per seed)
#   cf_*                  dev/simstudies/simstudy-carefulfit.R and
#                         simstudy-gaptol.R (same design; the latent process is
#                         simulated here by exact discretisation instead of
#                         ctGenerate(backend = 'r'), so the draws differ)
#   gA*, gB*, gC*         review/LAPLACE-gated-gaps-2026-09-24.md families A-C
#                         (the gaps job's defs2.R, copied, same draws)
#   gD*, gN*, anomS*      the same note's gap 3 (defs3.R, copied, same draws)
#   acnonlin              tests/testthat/test-julia-laplace-autocorrect.R
#   mvmix                 tests/testthat/test-julia-multivariate-mixed.R
#                         (ctGenerate replaced as for cf_*)
#   jflat                 tests/testthat/test-julia-convergence.R
#
# `bench_problem(model)` returns list(data = function(dataseed), model =
# function(), routes = allowed routes, idcols = id column(s)). Nothing in here touches the optimiser.

expmA <- function(A) as.matrix(Matrix::expm(A))

# ---- stochopt regimes (dev/stochopt/models.R) --------------------------------

sim_linear <- function(nsub, nobs, drift, diffchol, cintmean, cintsd, t0cov,
  dtrange = c(0.5, 1.5)) {
  k <- nrow(drift)
  Q <- diffchol %*% t(diffchol)
  Qinf <- matrix(solve(kronecker(diag(k), drift) + kronecker(drift, diag(k)),
    -as.vector(Q)), k)
  cache <- list()
  lapply(seq_len(nsub), function(i) {
    cint <- cintmean + rnorm(k, 0, cintsd)
    dts <- round(runif(nobs - 1, dtrange[1], dtrange[2]), 1)
    x <- matrix(NA, nobs, k)
    x[1, ] <- as.vector(t(chol(t0cov)) %*% rnorm(k))
    for (t in 2:nobs) {
      key <- as.character(dts[t - 1])
      if (is.null(cache[[key]])) {
        Ad <- expmA(drift * dts[t - 1])
        cache[[key]] <<- list(Ad = Ad,
          Lq = t(chol(Qinf - Ad %*% Qinf %*% t(Ad) + diag(1e-12, k))),
          B = solve(drift, Ad - diag(k)))
      }
      cc <- cache[[key]]
      x[t, ] <- as.vector(cc$Ad %*% x[t - 1, ] + cc$B %*% cint + cc$Lq %*% rnorm(k))
    }
    list(x = x, time = cumsum(c(0, dts)), id = i)
  })
}

measure <- function(sims, LAM, sd) {
  p <- nrow(LAM)
  do.call(rbind, lapply(sims, function(s) {
    y <- s$x %*% t(LAM) + matrix(rnorm(length(s$time) * p, 0, sd), ncol = p)
    colnames(y) <- paste0("Y", seq_len(p))
    data.frame(id = s$id, time = s$time, y)
  }))
}

diagchar <- function(names) {
  m <- matrix("0", length(names), length(names)); diag(m) <- names; m
}

invlog <- function(x) 1 / (1 + exp(-x))
drawcat <- function(eta, tau) {
  cum <- sapply(seq_along(tau), function(k) invlog(tau[k] - eta))
  p <- cbind(cum, 1)
  p <- cbind(p[, 1, drop = FALSE], t(apply(p, 1, diff)))
  apply(p, 1, function(pr) sample.int(length(pr), 1, prob = pmax(pr, 0)))
}

linear_model <- function(nlat, nman, LAMspec, re = TRUE) {
  lat <- paste0("eta", seq_len(nlat))
  m <- suppressMessages(ctModel(type = "ct", n.latent = nlat, n.manifest = nman,
    manifestNames = paste0("Y", seq_len(nman)), latentNames = lat,
    LAMBDA = LAMspec, MANIFESTMEANS = matrix(0, nman, 1),
    CINT = matrix(paste0("cint", seq_len(nlat))),
    MANIFESTVAR = diagchar(paste0("mv", seq_len(nman)))))
  m$pars$indvarying <- if (re) m$pars$param %in% paste0("cint", seq_len(nlat)) else FALSE
  m
}

stochopt_data <- function(model, seed) {
  set.seed(seed)
  if (model %in% c("panel", "panel5k", "small")) {
    N <- switch(model, panel = 1000, panel5k = 5000, small = 30)
    TT <- if (model == "small") 10 else 15
    drift <- matrix(c(-0.5, 0.2, -0.1, -0.4), 2)
    sims <- sim_linear(N, TT, drift, matrix(c(0.8, 0.2, 0, 0.6), 2),
      c(0.3, -0.2), c(0.5, 0.4), diag(1, 2))
    if (model != "small") {
      LAM <- matrix(c(1, 0.8, 0, 0, 0, 0, 1, 1.2), 4, 2)
      return(measure(sims, LAM, 0.5))
    }
    return(measure(sims, diag(2), 0.5))
  }
  if (model == "long") {
    drift <- matrix(c(-0.4, 0.15, 0, -0.1, -0.6, 0.2, 0.05, 0, -0.3), 3)
    sims <- sim_linear(1, 1500, drift, diag(c(0.7, 0.6, 0.5)), c(0.2, 0, -0.1),
      0, diag(1, 3), dtrange = c(0.2, 0.6))
    return(measure(sims, diag(3), 0.4))
  }
  if (model == "bigp") {
    k <- 4
    drift <- -diag(c(0.5, 0.4, 0.6, 0.3))
    drift[2, 1] <- 0.2; drift[3, 2] <- 0.15; drift[4, 3] <- -0.1; drift[1, 4] <- 0.1
    sims <- sim_linear(300, 20, drift, diag(0.6, k) + 0.1 * lower.tri(diag(k)),
      rep(0, k), rep(0.4, k), diag(1, k))
    LAM <- matrix(0, 8, k)
    for (j in 1:k) { LAM[2 * j - 1, j] <- 1; LAM[2 * j, j] <- 0.9 }
    return(measure(sims, LAM, 0.5))
  }
  if (model == "ordinal") {
    TAU <- c(-1.0, 0.4, 1.9)
    sims <- sim_linear(300, 10, matrix(-0.3), matrix(sqrt(0.8)), 0, 0.5, matrix(1))
    return(do.call(rbind, lapply(sims, function(s) {
      eta <- s$x[, 1]
      data.frame(id = s$id, time = s$time, o1 = drawcat(eta, TAU),
        o2 = drawcat(eta, TAU), o3 = drawcat(eta, TAU))
    })))
  }
  if (model == "nonlin") {
    # dx = (a x - b x^3 + c) dt + g dW : a double well at a > 0
    a <- 0.5; b <- 0.4; cc <- 0.1; g <- 0.6
    N <- 100; TT <- 30
    return(do.call(rbind, lapply(seq_len(N), function(i) {
      dts <- round(runif(TT - 1, 0.5, 1.5), 1)
      x <- numeric(TT); x[1] <- rnorm(1, 0, 1)
      for (t in 2:TT) {
        h <- dts[t - 1] / 50; z <- x[t - 1]
        for (s in 1:50) z <- z + (a * z - b * z^3 + cc) * h + g * sqrt(h) * rnorm(1)
        x[t] <- z
      }
      data.frame(id = i, time = cumsum(c(0, dts)),
        Y1 = x + rnorm(TT, 0, 0.3), Y2 = 0.8 * x + rnorm(TT, 0, 0.3))
    })))
  }
  stop("unknown stochopt model ", model)
}

stochopt_model <- function(model) {
  if (model %in% c("panel", "panel5k"))
    return(linear_model(2, 4, matrix(c(1, "l21", 0, 0, 0, 0, 1, "l42"), 4, 2)))
  if (model == "small") return(linear_model(2, 2, diag(2)))
  if (model == "long") return(linear_model(3, 3, diag(3), re = FALSE))
  if (model == "bigp") {
    k <- 4
    LAMspec <- matrix("0", 8, k)
    for (j in 1:k) { LAMspec[2 * j - 1, j] <- "1"; LAMspec[2 * j, j] <- paste0("l", j) }
    return(linear_model(k, 8, LAMspec))
  }
  if (model == "ordinal") {
    m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
      n.manifest = 3, manifestNames = c("o1", "o2", "o3"), latentNames = "eta1",
      LAMBDA = matrix(1, 3, 1), MANIFESTMEANS = matrix(0, 3, 1),
      CINT = matrix("cint"), T0MEANS = matrix(0), MANIFESTVAR = diag(0, 3),
      manifesttype = c(2L, 2L, 2L), ncategories = c(4L, 4L, 4L))))
    m$pars$indvarying <- m$pars$param %in% "cint"
    return(m)
  }
  if (model == "nonlin") {
    m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 2,
      manifestNames = c("Y1", "Y2"), latentNames = "eta1",
      LAMBDA = matrix(c(1, "l2"), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
      PARS = matrix(c("pa", "pb"), 1, 2),
      DRIFT = matrix("PARS[1,1] - PARS[1,2] * eta1^2"),
      CINT = matrix("cint"), DIFFUSION = matrix("diff"),
      MANIFESTVAR = diagchar(c("mv1", "mv2"))))
    m$pars$indvarying <- FALSE
    return(m)
  }
  stop("unknown stochopt model ", model)
}

# ---- carefulfit / gaptol cells (dev/simstudies) ------------------------------
#
# One latent, DRIFT -0.3, DIFFUSION 0.8 (Cholesky), T0VAR 1, a random CINT with
# sd 0.5, 60 subjects x 10 unit-spaced observations, and one of four measurement
# models. The studies drew the latent path with ctGenerate(backend = 'r'); it is
# the same process discretised exactly here.

CF_MEASURES <- list(
  gaussian = list(names = c("y1", "y2"), type = c(0L, 0L), ncat = c(0L, 0L)),
  binary   = list(names = c("b1", "b2", "b3"), type = c(1L, 1L, 1L),
    ncat = c(0L, 0L, 0L)),
  ordinal  = list(names = c("o1", "o2", "o3"), type = c(2L, 2L, 2L),
    ncat = c(4L, 4L, 4L)),
  mixed    = list(names = c("o1", "b1", "y1"), type = c(2L, 1L, 0L),
    ncat = c(4L, 0L, 0L)))

cf_data <- function(measure, seed, nsub = 60L, nobs = 10L) {
  set.seed(seed)
  spec <- CF_MEASURES[[measure]]
  drift <- -0.3; diffusion <- 0.8
  a <- exp(drift)
  b <- (a - 1) / drift
  q <- diffusion^2 * (a^2 - 1) / (2 * drift)
  cints <- stats::rnorm(nsub, 0, 0.5)
  d <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    eta <- numeric(nobs)
    eta[1] <- stats::rnorm(1, 0, 1)
    for (t in 2:nobs) eta[t] <- a * eta[t - 1] + b * cints[i] +
      stats::rnorm(1, 0, sqrt(q))
    data.frame(id = i, time = seq_len(nobs) - 1, eta = eta)
  }))
  eta <- d$eta
  TAU <- c(-1.0, 0.4, 1.9)
  for (j in seq_along(spec$names)) {
    d[[spec$names[j]]] <- switch(as.character(spec$type[j]),
      "0" = eta + stats::rnorm(length(eta), 0, 0.5),
      "1" = stats::rbinom(length(eta), 1, invlog(eta)),
      "2" = drawcat(eta, TAU))
  }
  d$eta <- NULL
  d
}

cf_model <- function(measure) {
  spec <- CF_MEASURES[[measure]]
  n <- length(spec$names)
  mvar <- diag(0, n)
  for (i in seq_len(n)) if (spec$type[i] == 0L) mvar[i, i] <- "mvar"
  args <- list(type = "ct", n.latent = 1, n.manifest = n,
    manifestNames = spec$names, latentNames = "eta1",
    LAMBDA = matrix(1, n, 1), MANIFESTMEANS = matrix(0, n, 1),
    CINT = matrix("cint"), T0MEANS = matrix(0), MANIFESTVAR = mvar,
    manifesttype = spec$type)
  if (any(spec$type == 2L)) args$ncategories <- spec$ncat
  m <- suppressWarnings(suppressMessages(do.call(ctModel, args)))
  m$pars$indvarying <- FALSE
  m$pars$indvarying[m$pars$param %in% "cint"] <- TRUE
  m
}

# ---- gated-gaps families A, B, C (defs2.R) -----------------------------------

gg_modelA <- function() suppressMessages(ctModel(silent = TRUE, type = "ct",
  CINT = "cint", MANIFESTMEANS = 0, LAMBDA = matrix(1),
  DRIFT = "drift|-log1p_exp(-param)|TRUE"))
gg_modelB <- function() suppressMessages(ctModel(silent = TRUE, type = "ct",
  CINT = "cint", MANIFESTMEANS = 0, LAMBDA = matrix(1),
  DRIFT = "drift|-log1p_exp(-param)",
  DIFFUSION = "diff|log1p_exp(param)|TRUE"))
gg_modelC <- function() suppressMessages(ctModel(silent = TRUE, type = "ct",
  CINT = "cint", MANIFESTMEANS = 0, LAMBDA = matrix(1),
  DRIFT = "drift|-log1p_exp(-param)",
  MANIFESTVAR = "mvar|log1p_exp(param)|TRUE"))

gg_genA <- function(seed, nsub, ntimes, mu, sdp, cor = TRUE) {
  set.seed(seed)
  baseline <- stats::rnorm(nsub, 2, 2)
  start <- stats::rnorm(nsub, baseline / 2, 1)
  raw <- if (cor) stats::rnorm(nsub, mu + (baseline - 2) / 2, sdp) else stats::rnorm(nsub, mu, sdp)
  drift <- -log1p(exp(-raw))
  rows <- lapply(seq_len(nsub), function(i) {
    a <- drift[i]; decay <- exp(a)
    intercept <- (baseline[i] / a) * (decay - 1)
    innovation <- sqrt(0.25 * (exp(2 * a) - 1) / (2 * a))
    latent <- numeric(ntimes); latent[1] <- start[i]
    for (t in seq_len(ntimes - 1L)) latent[t + 1L] <- decay * latent[t] + intercept + stats::rnorm(1, 0, innovation)
    data.frame(id = i, time = seq_len(ntimes) - 1L, Y1 = latent + stats::rnorm(ntimes, 0, 0.5))
  })
  do.call(rbind, rows)
}

gg_genBC <- function(family, seed, nsub, ntimes, mu, sdp) {
  set.seed(seed)
  a <- -0.5; decay <- exp(a)
  baseline <- stats::rnorm(nsub, 2, 1.5)
  start <- stats::rnorm(nsub, 2, 1)
  raw <- stats::rnorm(nsub, mu, sdp)
  sdv <- log1p(exp(raw))
  rows <- lapply(seq_len(nsub), function(i) {
    diffsd <- if (family == "B") sdv[i] else 0.6
    noisesd <- if (family == "C") sdv[i] else 0.4
    innovation <- sqrt(diffsd^2 * (exp(2 * a) - 1) / (2 * a))
    latent <- numeric(ntimes); latent[1] <- start[i]
    for (t in seq_len(ntimes - 1L)) latent[t + 1L] <- decay * latent[t] +
      baseline[i] * (1 - decay) + stats::rnorm(1, 0, innovation)
    data.frame(id = i, time = seq_len(ntimes) - 1L, Y1 = latent + stats::rnorm(ntimes, 0, noisesd))
  })
  do.call(rbind, rows)
}

# Family A's configs were kept in the gaps job's thr2.rds; this is the same
# grid, in the same order (checked against that file).
gg_configsA <- expand.grid(seed = 1:2, ntimes = c(6L, 12L), mu = c(1, 3),
  sdp = c(1, 2.5))
gg_configsA$nsub <- 40L
gg_configsBC <- expand.grid(seed = 1:2, ntimes = c(6L, 12L), mu = c(0, -1),
  stringsAsFactors = FALSE)
gg_configsBC$sdp <- ifelse(gg_configsBC$mu == 0, 1, 1.5)
gg_configsBC$nsub <- 40L

# ---- gated-gaps gap 3: N (nested), D (binary), S (AnomAuth) (defs3.R) --------

gg_modelN <- function() {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", manifestNames = "Y1",
    latentNames = "eta1", LAMBDA = matrix(1), MANIFESTMEANS = matrix(0),
    CINT = matrix("cint"), DRIFT = matrix("drift|-log1p_exp(-param)"),
    T0VAR = matrix(0.5), id = c("subject", "study"), silent = TRUE)))
  m$pars$indvarying <- m$pars$param %in% c("drift", "cint")
  m$pars$indvarying_study <- m$pars$param %in% "cint"
  m
}

gg_genN <- function(seed, nstudy, npersub, ntimes, mu, sdp) {
  set.seed(seed)
  rows <- list(); sid <- 0
  for (g in seq_len(nstudy)) {
    studycint <- stats::rnorm(1, 0, 1)
    for (j in seq_len(npersub)) {
      sid <- sid + 1
      raw <- stats::rnorm(1, mu, sdp)
      a <- -log1p(exp(-raw)); decay <- exp(a)
      baseline <- 2 + studycint + stats::rnorm(1, 0, 1)
      intercept <- (baseline / a) * (decay - 1)
      innovation <- sqrt(0.25 * (exp(2 * a) - 1) / (2 * a))
      latent <- numeric(ntimes); latent[1] <- stats::rnorm(1, 0, 0.7)
      for (t in seq_len(ntimes - 1L)) latent[t + 1L] <- decay * latent[t] + intercept +
        stats::rnorm(1, 0, innovation)
      rows[[length(rows) + 1L]] <- data.frame(subject = sid, study = g,
        time = seq_len(ntimes) - 1L, Y1 = latent + stats::rnorm(ntimes, 0, 0.5))
    }
  }
  do.call(rbind, rows)
}
gg_configsN <- data.frame(seed = c(1, 2, 1, 2), nstudy = 8L, npersub = 5L,
  ntimes = c(6L, 6L, 10L, 10L), mu = 1, sdp = 1.5)

gg_modelD <- function() {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = c("B1", "B2", "B3"), latentNames = "eta1",
    LAMBDA = matrix(1, 3, 1), manifesttype = c(1, 1, 1),
    MANIFESTMEANS = matrix(c(0, "thr2", "thr3"), 3, 1), MANIFESTVAR = matrix(0, 3, 3),
    CINT = matrix("cint||TRUE"), DRIFT = matrix("drift|-log1p_exp(-param)|TRUE"),
    T0VAR = matrix(0.5), silent = TRUE)))
  m$pars$indvarying <- m$pars$param %in% c("drift", "cint")
  m
}

gg_genD <- function(seed, nsub, ntimes, mu, sdp) {
  set.seed(seed)
  rows <- lapply(seq_len(nsub), function(i) {
    raw <- stats::rnorm(1, mu, sdp); a <- -log1p(exp(-raw)); decay <- exp(a)
    baseline <- stats::rnorm(1, 0, 1)
    intercept <- (baseline / a) * (decay - 1)
    innovation <- sqrt(1 * (exp(2 * a) - 1) / (2 * a))
    latent <- numeric(ntimes); latent[1] <- stats::rnorm(1, 0, 0.7)
    for (t in seq_len(ntimes - 1L)) latent[t + 1L] <- decay * latent[t] + intercept +
      stats::rnorm(1, 0, innovation)
    th <- c(0, 0.5, -0.5)
    Y <- sapply(th, function(m) stats::rbinom(ntimes, 1, stats::plogis(latent + m)))
    data.frame(id = i, time = seq_len(ntimes) - 1L, B1 = Y[, 1], B2 = Y[, 2], B3 = Y[, 3])
  })
  do.call(rbind, rows)
}
gg_configsD <- data.frame(seed = c(1, 2, 1, 2), nsub = 60L, ntimes = c(8L, 8L, 16L, 16L),
  mu = 1, sdp = 1.5)

# AnomAuth (shipped with ctsem): anomia alone, subjects with at least 3 observed
# waves, 800 of them at random; random CINT and a random -log1p_exp auto-effect.
gg_modelS <- function() {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", manifestNames = "Y1",
    latentNames = "anom", LAMBDA = matrix(1), MANIFESTMEANS = matrix(0),
    CINT = matrix("cint||TRUE"), DRIFT = matrix("drift|-log1p_exp(-param)|TRUE"),
    T0MEANS = matrix("t0m"), silent = TRUE)))
  m$pars$indvarying <- m$pars$param %in% c("drift", "cint")
  m
}
gg_genS <- function(seed, nsub) {
  e <- new.env(); utils::data("AnomAuth", package = "ctsem", envir = e)
  long <- suppressMessages(ctWideToLong(e$AnomAuth, Tpoints = 5, n.manifest = 2,
    manifestNames = c("Y1", "Y2")))
  long <- as.data.frame(suppressMessages(ctDeintervalise(long)))
  long <- long[!is.na(long$Y1), c("id", "time", "Y1")]
  nobs <- table(long$id); keep <- as.numeric(names(nobs)[nobs >= 3])
  set.seed(seed)
  ids <- sort(sample(keep, min(nsub, length(keep))))
  long[long$id %in% ids, ]
}
gg_configsS <- data.frame(seed = c(1, 2), nsub = 800L)

# ---- test fixtures -------------------------------------------------------------

# test-julia-laplace-autocorrect.R: one random effect on a -log1p_exp DRIFT.
ac_model <- function() {
  m <- suppressMessages(ctModel(silent = TRUE, type = "ct", CINT = 0,
    MANIFESTMEANS = 0, LAMBDA = matrix(1), T0MEANS = matrix(0),
    DRIFT = "drift|-log1p_exp(-param)|TRUE"))
  m$pars$indvarying <- m$pars$param %in% "drift"
  m
}
ac_data <- function(seed = 3L, nsubjects = 40L, ntimes = 8L) {
  set.seed(seed)
  drift <- -log1p(exp(-stats::rnorm(nsubjects, 1, 1)))
  rows <- lapply(seq_len(nsubjects), function(i) {
    a <- drift[i]; decay <- exp(a)
    innovation <- sqrt(0.25 * (exp(2 * a) - 1) / (2 * a))
    latent <- numeric(ntimes); latent[1] <- stats::rnorm(1, 0, 1)
    for (t in seq_len(ntimes - 1L)) latent[t + 1L] <- decay * latent[t] +
      stats::rnorm(1, 0, innovation)
    data.frame(id = i, time = seq_len(ntimes) - 1L,
      Y1 = latent + stats::rnorm(ntimes, 0, 0.3))
  })
  do.call(rbind, rows)
}

# test-julia-multivariate-mixed.R: two latents, ordinal + binary indicators on
# the first and Gaussian on the second, random CINTs on both -- a random-effect
# covariance the data cannot fully separate (the "rank-deficient ridge"). The
# test draws the latent path with ctGenerate(backend = 'r'); here it is the
# same process, discretised exactly at unit intervals.
mvmix_data <- function(seed = 404L, nsubjects = 30L, nobs = 8L) {
  set.seed(seed)
  tau <- c(-0.8, 0.3, 1.4)
  drawcat3 <- function(eta) drawcat(eta, tau)
  drift <- matrix(c(-0.4, 0.15, -0.10, -0.6), 2, 2, byrow = TRUE)
  dchol <- matrix(c(0.7, 0, 0.2, 0.6), 2, 2)
  Q <- dchol %*% t(dchol)
  Ad <- expmA(drift)
  Qinf <- matrix(solve(kronecker(diag(2), drift) + kronecker(drift, diag(2)),
    -as.vector(Q)), 2)
  Lq <- t(chol(Qinf - Ad %*% Qinf %*% t(Ad)))
  B <- solve(drift, Ad - diag(2))
  cint1 <- stats::rnorm(nsubjects, 0, 0.4)
  cint2 <- stats::rnorm(nsubjects, 0, 0.4)
  d <- do.call(rbind, lapply(seq_len(nsubjects), function(i) {
    x <- matrix(NA_real_, nobs, 2)
    x[1, ] <- stats::rnorm(2)
    cint <- c(cint1[i], cint2[i])
    for (t in 2:nobs) x[t, ] <- as.vector(Ad %*% x[t - 1, ] + B %*% cint +
      Lq %*% stats::rnorm(2))
    data.frame(id = i, time = seq_len(nobs) - 1, e1 = x[, 1], e2 = x[, 2])
  }))
  d$o1 <- drawcat3(d$e1)
  d$o2 <- drawcat3(d$e1)
  d$b1 <- stats::rbinom(nrow(d), 1, invlog(d$e1))
  d$y1 <- d$e2 + stats::rnorm(nrow(d), 0, 0.4)
  d$y2 <- d$e2 + stats::rnorm(nrow(d), 0, 0.4)
  d$e1 <- NULL
  d$e2 <- NULL
  d
}
mvmix_model <- function() {
  lambda <- matrix(0, 5, 2)
  lambda[1:3, 1] <- 1
  lambda[4:5, 2] <- 1
  mvar <- diag(0, 5)
  mvar[4, 4] <- "mv1"
  mvar[5, 5] <- "mv2"
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 2,
    n.manifest = 5, manifestNames = c("o1", "o2", "b1", "y1", "y2"),
    latentNames = c("eta1", "eta2"), manifesttype = c(2L, 2L, 1L, 0L, 0L),
    ncategories = c(4L, 4L, 0L, 0L, 0L), LAMBDA = lambda,
    MANIFESTMEANS = matrix(0, 5, 1),
    CINT = matrix(c("cint1", "cint2"), 2, 1), T0MEANS = matrix(0, 2, 1),
    MANIFESTVAR = mvar)))
  m$pars$indvarying <- FALSE
  m$pars$indvarying[m$pars$param %in% c("cint1", "cint2")] <- TRUE
  m
}

# test-julia-convergence.R: the flat-transform fixture (25 x 25, seed 3), whose
# point is the stored start with drift at raw 8 (starts.R, `flatdrift8`).
jflat_data <- function(seed = 3L, nsubjects = 25L, ntimes = 25L) {
  set.seed(seed)
  baseline <- stats::rnorm(nsubjects, 2, 2)
  start <- stats::rnorm(nsubjects, baseline / 2, 1)
  drift <- -log1p(exp(-stats::rnorm(nsubjects, baseline / 2, 0.5)))
  diffusion <- 0.5
  noise <- 0.5
  rows <- lapply(seq_len(nsubjects), function(i) {
    a <- drift[i]
    decay <- exp(a)
    intercept <- (baseline[i] / a) * (decay - 1)
    innovation <- sqrt(diffusion^2 * (exp(2 * a) - 1) / (2 * a))
    latent <- numeric(ntimes)
    latent[1] <- start[i]
    for (t in seq_len(ntimes - 1L)) {
      latent[t + 1L] <- decay * latent[t] + intercept +
        stats::rnorm(1, 0, innovation)
    }
    data.frame(id = i, time = seq_len(ntimes) - 1L,
      Y1 = latent + stats::rnorm(ntimes, 0, noise),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
jflat_model <- function() suppressMessages(ctModel(silent = TRUE,
  type = "ct", CINT = "cint", MANIFESTMEANS = 0, LAMBDA = matrix(1),
  DRIFT = "drift|-log1p_exp(-param)|TRUE"))

# ---- the registry ----------------------------------------------------------------

bench_problem <- function(model) {
  fixed <- function(gen) function(dataseed) gen()  # config fixes its own data
  mk <- function(data, build, routes, idcol = "id", reference = TRUE,
    datafixed = FALSE, note = "") list(data = data, model = build,
    routes = routes, idcol = idcol, reference = reference,
    datafixed = datafixed, note = note)
  if (model %in% c("panel", "panel5k", "long", "ordinal", "nonlin", "small", "bigp")) {
    route <- if (model == "ordinal") "laplace" else "augmented"
    return(mk(function(s) stochopt_data(model, as.integer(s)),
      function() stochopt_model(model), c(route, "augmented", "laplace", "auto")))
  }
  if (grepl("^cf_", model)) {
    measure <- sub("^cf_", "", model)
    if (!measure %in% names(CF_MEASURES)) stop("unknown cf measure ", measure)
    return(mk(function(s) cf_data(measure, as.integer(s)),
      function() cf_model(measure), c("augmented", "laplace", "auto")))
  }
  if (grepl("^g[ABC][0-9]+$", model)) {
    fam <- substr(model, 2, 2); k <- as.integer(substring(model, 3))
    if (fam == "A") {
      cf <- gg_configsA[k, ]
      if (anyNA(cf$seed)) stop("no config ", model)
      return(mk(fixed(function() gg_genA(cf$seed, cf$nsub, cf$ntimes, cf$mu, cf$sdp)),
        gg_modelA, c("laplace", "augmented", "auto"), datafixed = TRUE))
    }
    cf <- gg_configsBC[k, ]
    if (anyNA(cf$seed)) stop("no config ", model)
    return(mk(fixed(function() gg_genBC(fam, cf$seed, cf$nsub, cf$ntimes, cf$mu, cf$sdp)),
      if (fam == "B") gg_modelB else gg_modelC, c("laplace", "augmented", "auto"),
      datafixed = TRUE))
  }
  if (grepl("^gD[0-9]+$", model)) {
    cf <- gg_configsD[as.integer(substring(model, 3)), ]
    if (anyNA(cf$seed)) stop("no config ", model)
    return(mk(fixed(function() gg_genD(cf$seed, cf$nsub, cf$ntimes, cf$mu, cf$sdp)),
      gg_modelD, c("laplace", "augmented", "auto"), datafixed = TRUE))
  }
  if (grepl("^gN[0-9]+$", model)) {
    cf <- gg_configsN[as.integer(substring(model, 3)), ]
    if (anyNA(cf$seed)) stop("no config ", model)
    return(mk(fixed(function() gg_genN(cf$seed, cf$nstudy, cf$npersub, cf$ntimes,
      cf$mu, cf$sdp)), gg_modelN, c("laplace", "auto"), idcol = "study",
      datafixed = TRUE))
  }
  if (grepl("^anomS[12]$", model)) {
    cf <- gg_configsS[as.integer(substring(model, 6)), ]
    return(mk(fixed(function() gg_genS(cf$seed, cf$nsub)), gg_modelS,
      c("laplace", "augmented", "auto"), datafixed = TRUE))
  }
  if (model == "acnonlin") return(mk(fixed(function() ac_data()), ac_model,
    c("laplace", "augmented", "auto"), datafixed = TRUE))
  if (model == "mvmix") return(mk(fixed(function() mvmix_data()), mvmix_model,
    c("laplace", "augmented", "auto"), datafixed = TRUE))
  if (model == "jflat") return(mk(fixed(function() jflat_data()), jflat_model,
    c("laplace", "augmented", "auto"), datafixed = TRUE))
  stop("unknown bench model '", model, "'")
}

# Stored once per (model, data) in `store`, written atomically, and read back
# from there by every later cell and every build.
bench_data <- function(model, dataseed, store) {
  P <- bench_problem(model)
  key <- if (isTRUE(P$datafixed)) "cfg" else as.character(dataseed)
  dir.create(store, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(store, sprintf("%s-%s.rds", model, key))
  if (!file.exists(f)) {
    d <- P$data(dataseed)
    tmp <- paste0(f, ".tmp", Sys.getpid())
    saveRDS(d, tmp)
    file.rename(tmp, f)
  }
  list(data = readRDS(f), file = f, md5 = unname(tools::md5sum(f)))
}
