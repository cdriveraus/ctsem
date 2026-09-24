# Regime-balanced model set for the optimiser prototypes. Own simulators with
# exact discretisation (Euler substeps for the nonlinear one), so nothing here
# depends on ctGenerate's draw stream.
#
#   panel    2 latents, 4 indicators, random CINTs, N=1000, T=15   augmented
#   panel5k  the same with N=5000 (demonstration of batching and finish)
#   long     3 latents, 3 indicators, N=1, T=1500, full DRIFT         filter
#   ordinal  1 latent, 3 ordinal indicators, random CINT, N=300, T=10 laplace
#   nonlin   1 latent cubic drift, 2 indicators, N=100, T=30           EKF
#   small    2 latents, 2 indicators, random CINTs, N=30, T=10       augmented
#   bigp     4 latents, 8 indicators, full DRIFT/DIFFUSION, random CINTs,
#            N=300, T=20                                             augmented

expmA <- function(A) as.matrix(Matrix::expm(A))

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

linear_model <- function(nlat, nman, LAMspec, re = TRUE, fulldrift = TRUE) {
  lat <- paste0("eta", seq_len(nlat))
  m <- suppressMessages(ctModel(type = "ct", n.latent = nlat, n.manifest = nman,
    manifestNames = paste0("Y", seq_len(nman)), latentNames = lat,
    LAMBDA = LAMspec, MANIFESTMEANS = matrix(0, nman, 1),
    CINT = matrix(paste0("cint", seq_len(nlat))),
    MANIFESTVAR = diagchar(paste0("mv", seq_len(nman)))))
  m$pars$indvarying <- if (re) m$pars$param %in% paste0("cint", seq_len(nlat)) else FALSE
  m
}

make_problem <- function(model, seed) {
  set.seed(seed)
  if (model %in% c("panel", "panel5k", "small")) {
    N <- switch(model, panel = 1000, panel5k = 5000, small = 30)
    TT <- if (model == "small") 10 else 15
    drift <- matrix(c(-0.5, 0.2, -0.1, -0.4), 2)
    sims <- sim_linear(N, TT, drift, matrix(c(0.8, 0.2, 0, 0.6), 2),
      c(0.3, -0.2), c(0.5, 0.4), diag(1, 2))
    if (model != "small") {
      LAM <- matrix(c(1, 0.8, 0, 0, 0, 0, 1, 1.2), 4, 2)
      d <- measure(sims, LAM, 0.5)
      m <- linear_model(2, 4, matrix(c(1, "l21", 0, 0, 0, 0, 1, "l42"), 4, 2))
    } else {
      d <- measure(sims, diag(2), 0.5)
      m <- linear_model(2, 2, diag(2))
    }
    return(list(d = d, m = m, route = "augmented"))
  }
  if (model == "long") {
    drift <- matrix(c(-0.4, 0.15, 0, -0.1, -0.6, 0.2, 0.05, 0, -0.3), 3)
    sims <- sim_linear(1, 1500, drift, diag(c(0.7, 0.6, 0.5)), c(0.2, 0, -0.1),
      0, diag(1, 3), dtrange = c(0.2, 0.6))
    d <- measure(sims, diag(3), 0.4)
    m <- linear_model(3, 3, diag(3), re = FALSE)
    return(list(d = d, m = m, route = "augmented"))
  }
  if (model == "bigp") {
    k <- 4
    drift <- -diag(c(0.5, 0.4, 0.6, 0.3))
    drift[2, 1] <- 0.2; drift[3, 2] <- 0.15; drift[4, 3] <- -0.1; drift[1, 4] <- 0.1
    sims <- sim_linear(300, 20, drift, diag(0.6, k) + 0.1 * lower.tri(diag(k)),
      rep(0, k), rep(0.4, k), diag(1, k))
    LAM <- matrix(0, 8, k)
    for (j in 1:k) { LAM[2 * j - 1, j] <- 1; LAM[2 * j, j] <- 0.9 }
    d <- measure(sims, LAM, 0.5)
    LAMspec <- matrix("0", 8, k)
    for (j in 1:k) { LAMspec[2 * j - 1, j] <- "1"; LAMspec[2 * j, j] <- paste0("l", j) }
    m <- linear_model(k, 8, LAMspec)
    return(list(d = d, m = m, route = "augmented"))
  }
  if (model == "ordinal") {
    TAU <- c(-1.0, 0.4, 1.9)
    sims <- sim_linear(300, 10, matrix(-0.3), matrix(sqrt(0.8)), 0, 0.5, matrix(1))
    d <- do.call(rbind, lapply(sims, function(s) {
      eta <- s$x[, 1]
      data.frame(id = s$id, time = s$time, o1 = drawcat(eta, TAU),
        o2 = drawcat(eta, TAU), o3 = drawcat(eta, TAU))
    }))
    m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
      n.manifest = 3, manifestNames = c("o1", "o2", "o3"), latentNames = "eta1",
      LAMBDA = matrix(1, 3, 1), MANIFESTMEANS = matrix(0, 3, 1),
      CINT = matrix("cint"), T0MEANS = matrix(0), MANIFESTVAR = diag(0, 3),
      manifesttype = c(2L, 2L, 2L), ncategories = c(4L, 4L, 4L))))
    m$pars$indvarying <- m$pars$param %in% "cint"
    return(list(d = d, m = m, route = "laplace"))
  }
  if (model == "nonlin") {
    # dx = (a x - b x^3 + c) dt + g dW : a double well at a > 0
    a <- 0.5; b <- 0.4; cc <- 0.1; g <- 0.6
    N <- 100; TT <- 30
    d <- do.call(rbind, lapply(seq_len(N), function(i) {
      dts <- round(runif(TT - 1, 0.5, 1.5), 1)
      x <- numeric(TT); x[1] <- rnorm(1, 0, 1)
      for (t in 2:TT) {
        h <- dts[t - 1] / 50; z <- x[t - 1]
        for (s in 1:50) z <- z + (a * z - b * z^3 + cc) * h + g * sqrt(h) * rnorm(1)
        x[t] <- z
      }
      data.frame(id = i, time = cumsum(c(0, dts)),
        Y1 = x + rnorm(TT, 0, 0.3), Y2 = 0.8 * x + rnorm(TT, 0, 0.3))
    }))
    m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 2,
      manifestNames = c("Y1", "Y2"), latentNames = "eta1",
      LAMBDA = matrix(c(1, "l2"), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
      PARS = matrix(c("pa", "pb"), 1, 2),
      DRIFT = matrix("PARS[1,1] - PARS[1,2] * eta1^2"),
      CINT = matrix("cint"), DIFFUSION = matrix("diff"),
      MANIFESTVAR = diagchar(c("mv1", "mv2"))))
    m$pars$indvarying <- FALSE
    return(list(d = d, m = m, route = "augmented"))
  }
  stop("unknown model ", model)
}
