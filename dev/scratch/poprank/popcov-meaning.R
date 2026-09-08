## What does a fixed POPCOV entry actually mean?
##
## The code says the diagonal is a natural-scale sd (converted by the
## transform's slope) and the off-diagonal "a correlation ... must lie in
## [-1,1]". But the free branch gives the correlation cell the transform
## 2/(1+exp(-param))-1, and the whole T0VAR then goes through sdcovsqrt2cov,
## which applies constraincorsqrt1 -- a row-normalising map. So the number in
## the cell is an *input* to that map, not the correlation itself. Algebra for
## the 2x2 case says a cell of 0.5 comes out near 0.79, and that a cell of 0
## comes out exactly 0. Checked here rather than asserted.
##
## poprank is off throughout: this is about POPCOV alone.
Sys.setenv(NOT_CRAN = 'true')
devtools::load_all('.', compile = FALSE, quiet = TRUE)
say <- function(...) { cat('@@', ..., '\n', sep = ''); flush.console() }
say('LOADED')

## two mean-affecting random effects, so nothing is partially identified
set.seed(4); nsub <- 120L; nt <- 6L
o <- vector('list', nsub)
for (i in seq_len(nsub)) {
  m1 <- rnorm(1, 0, .5); m2 <- rnorm(1, 0, .5)
  eta <- numeric(nt); eta[1] <- rnorm(1)
  for (t in 2:nt) eta[t] <- exp(-.5) * eta[t - 1] + rnorm(1, 0, .4)
  o[[i]] <- data.frame(id = i, time = seq_len(nt) - 1,
    Y1 = eta + m1 + rnorm(nt, 0, .3), Y2 = eta + m2 + rnorm(nt, 0, .3))
}
dat <- do.call(rbind, o)

mk <- function() {
  m <- suppressMessages(ctModel(type = 'ct', n.latent = 1, n.manifest = 2,
    LAMBDA = matrix(c(1, 1), 2, 1), manifestNames = c('Y1', 'Y2'),
    latentNames = 'eta', DRIFT = matrix(-0.5), DIFFUSION = matrix(0.4),
    CINT = matrix(0), T0MEANS = matrix(0), T0VAR = matrix(1),
    MANIFESTMEANS = matrix(c('mm1|param', 'mm2|param'), 2, 1),
    MANIFESTVAR = matrix(c('mv1|log1p_exp(param)', 0, 0,
      'mv2|log1p_exp(param)'), 2, 2)))
  m$pars$indvarying <- FALSE
  m$pars$indvarying[m$pars$matrix %in% 'MANIFESTMEANS'] <- TRUE
  m
}
prep <- function(m) ctsem:::ctModelStatesAndPARS(
  ctsem:::ctModel0DRIFT(m, TRUE)$pars, statenames = 'eta', tdprednames = NULL)

for (cell in c(NA, 0, 0.3, 0.5, 0.8)) {
  m <- mk()
  if (!is.na(cell)) {
    p <- prep(m); m$pars <- p
    m[['POPCOV']] <- ctsem:::.ctModelPopCov(p)
    nms <- rownames(m[['POPCOV']])
    m[['POPCOV']][nms[2], nms[1]] <- cell
  }
  set.seed(77)
  f <- try(suppressWarnings(suppressMessages(ctFit(datalong = dat, model = m,
    backend = 'julia', intoverpop = 'augmented', poprank = NA, cores = 1L,
    verbose = 0L, optimcontrol = list(estonly = TRUE)))), silent = TRUE)
  if (inherits(f, 'try-error')) {
    say(sprintf('POPCOV cell %-4s FIT FAILED: %s', cell,
      substr(conditionMessage(attr(f, 'condition')), 1, 90))); next }
  cv <- try(ctsem:::.ctBackendRawPopCov(f), silent = TRUE)
  if (inherits(cv, 'try-error') || is.null(cv)) {
    say(sprintf('POPCOV cell %-4s no covariance read', cell)); next }
  cov <- cv[[1]]$cov
  corr <- cov[2, 1] / sqrt(cov[1, 1] * cov[2, 2])
  say(sprintf('POPCOV cell %-4s -> implied population correlation %+.4f  (sds %.4f, %.4f)',
    if (is.na(cell)) 'free' else format(cell), corr, sqrt(cov[1, 1]), sqrt(cov[2, 2])))
}

say('')
say('and the diagonal, which the code converts by the transform slope:')
for (sd in c(NA, 0.3, 0.8)) {
  m <- mk()
  if (!is.na(sd)) {
    p <- prep(m); m$pars <- p
    m[['POPCOV']] <- ctsem:::.ctModelPopCov(p)
    nms <- rownames(m[['POPCOV']])
    m[['POPCOV']][nms[1], nms[1]] <- sd
  }
  set.seed(77)
  f <- try(suppressWarnings(suppressMessages(ctFit(datalong = dat, model = m,
    backend = 'julia', intoverpop = 'augmented', poprank = NA, cores = 1L,
    verbose = 0L, optimcontrol = list(estonly = TRUE)))), silent = TRUE)
  if (inherits(f, 'try-error')) { say('  sd ', sd, ' FAILED'); next }
  cv <- ctsem:::.ctBackendRawPopCov(f)
  say(sprintf('  POPCOV diagonal %-4s -> raw-scale sd %.4f',
    if (is.na(sd)) 'free' else format(sd), sqrt(cv[[1]]$cov[1, 1])))
}
say('DONE')
