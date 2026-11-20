# A reduced-rank level in blocks, as the hc-wellbeing reanalysis fits them: a
# subject's level in one block, item intercepts in another named by matrix,
# one of them about a stated mean. The intercepts are generated from two
# dimensions, which their block's one does not span, independently of the
# level. Shared by the score tests and the search that acts on them
# (test-ctScoreTest.R, test-ctModelImprove.R).
.score_blocks_fit <- function() fit_cached('scoretest blocks', {
  set.seed(21)
  lam <- c(1, 0.8, 1.2, 1)
  d <- do.call(rbind, lapply(seq_len(150L), function(i) {
    times <- cumsum(c(0, stats::runif(9L, 0.5, 1.5)))
    lev <- stats::rnorm(1, 0, 0.6)
    f <- stats::rnorm(2)
    mm <- c(0, 0.5 * f[1], 0.5 * f[1] + 0.5 * f[2], -0.5 * f[2])
    eta <- numeric(10L); eta[1] <- stats::rnorm(1)
    for (t in 2:10) {
      e <- exp(-0.6 * (times[t] - times[t - 1]))
      eta[t] <- e * eta[t - 1] + (1 - e) * lev / 0.6 +
        stats::rnorm(1, 0, sqrt((1 - e^2) / 1.2))
    }
    y <- sapply(1:4, function(j) lam[j] * eta + mm[j] + stats::rnorm(10L, 0, 0.4))
    data.frame(id = i, time = times, Y1 = y[, 1], Y2 = y[, 2], Y3 = y[, 3], Y4 = y[, 4])
  }))
  m <- suppressMessages(ctModel(type = 'ct', manifestNames = paste0('Y', 1:4),
    latentNames = 'eta', LAMBDA = matrix(lam), DRIFT = matrix('dr'),
    CINT = matrix(ctPar('lev', indvarying = TRUE)), DIFFUSION = matrix('q'),
    MANIFESTMEANS = matrix(c(ctPar('mm1', indvarying = TRUE, mean = 0),
      ctPar('mm2', indvarying = TRUE), ctPar('mm3', indvarying = TRUE),
      ctPar('mm4', indvarying = TRUE))), T0MEANS = matrix(0), T0VAR = matrix(1)))
  blocks <- list(dynamics = list(rank = 1, effects = 'lev'),
    items = list(rank = 1, effects = 'MANIFESTMEANS'))
  list(data = d, blocks = blocks, fit = suppressWarnings(suppressMessages(ctFit(d, m,
    backend = 'julia', intoverpop = 'laplace', popblocks = blocks, cores = 1,
    verbose = 0))))
})
